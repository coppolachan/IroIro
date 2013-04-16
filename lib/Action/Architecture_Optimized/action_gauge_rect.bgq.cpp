/*!
  @file action_gauge_rect.bgq.cpp
  @brief Definition of the ActionGaugeRect class

  Optimized version for BlueGeneQ architecture

  Time-stamp: <2013-04-16 15:41:25 cossu>
*/

#include "Action/action_gauge_rect.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/fieldUtils.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include "lib/Main/Geometry/mapping.hpp"
#include "include/messages_macros.hpp"

#include <omp.h>
#include "bgqthread.h"

using namespace SUNmatUtils;
using namespace FieldUtils;

typedef struct GaugeConfigType{
  __complex__ double v[NC_*NC_];
}GaugePtr;

double ActionGaugeRect::calc_H(){
  using namespace Mapping;

  const Staples stpl_;
  GaugeField1D Cup1, Cup2;
  GaugeField1D res;
  GaugeField1D UpNu, UpMu;
  GaugeField1D U_mu, U_nu;

  double plaqF = 0.0;
  double rectF = 0.0;

  //pointers
  double* UpMu_ptr = UpMu.data.getaddr(0);
  double* UpNu_ptr = UpNu.data.getaddr(0);
  double* res_ptr = res.data.getaddr(0);
  double* U_ptr = u_->data.getaddr(0);
  double* U_nu_ptr = u_->data.getaddr(0);
  double* U_mu_ptr = u_->data.getaddr(0);
  double* Cup1_ptr = Cup1.data.getaddr(0);
  double* Cup2_ptr = Cup2.data.getaddr(0);

  BGQThread_Init();

  // From Matsufuru-san code
#pragma omp parallel
  {
    const int nid = omp_get_num_threads();
    const int tid = omp_get_thread_num();
    const int is = tid*Nvol_/nid;
    const int ie = (tid + 1)*Nvol_/nid;
    const int ns = ie - is;
    const int CC = NC_*NC_;
    const int str = is*CC;
    const int CC2 = 2*NC_*NC_;
    const int str2 = is*CC2;
    
    for(int mu=0; mu<NDIM_; ++mu){
      U_mu_ptr = U_ptr+mu*Nvol_*CC2; // assumes contiguity of U_mu     
      for(int nu = mu+1; nu<NDIM_; ++nu){
	U_nu_ptr = U_ptr+nu*Nvol_*CC2;      

	shiftField(UpNu,U_mu_ptr,nu,Forward());
	shiftField(UpMu,U_nu_ptr,mu,Forward());
	
	BGWilsonSU3_MatMult_NND(Cup1_ptr+str2,U_nu_ptr+str2,UpNu_ptr+str2,
				UpMu_ptr+str2,ns);
	BGWilsonSU3_MatMult_NND(Cup2_ptr+str2,U_mu_ptr+str2,UpMu_ptr+str2,
				UpNu_ptr+str2,ns);

	// plaquette term
	BGWilsonSU3_MatMult_ND(res_ptr+str2,U_mu_ptr+str2,Cup1_ptr+str2,ns);

#pragma omp for reduction(+:plaqF) //nb NOT "parallel for" just "for"
	for(int site=0; site<Nvol_; ++site)
	  plaqF += res_ptr[CC2*site]+res_ptr[8+CC2*site]+res_ptr[16+CC2*site]; ///ReTr

	
	shiftField(UpMu, Cup2_ptr,mu,Forward());      //Cup2(x+mu) 

	BGWilsonSU3_MatMult_NND(res_ptr+str2,
				U_nu_ptr+str2,UpNu_ptr+str2,UpMu_ptr+str2,ns);
	BGWilsonSU3_MatMult_ND(UpNu_ptr+str2,U_mu_ptr+str2,res_ptr+str2,ns);

#pragma omp for reduction(+:rectF)
	for(int site=0; site<Nvol_; ++site){
	  rectF += UpNu_ptr[CC2*site]+ UpNu_ptr[8+CC2*site]+ UpNu_ptr[16+CC2*site];///ReTr
	}
	
	shiftField(UpNu,Cup1_ptr,nu,Forward());   //Cup1(x+nu)
	shiftField(UpMu,U_nu_ptr,mu,Forward());   //U_nu(x+mu)

	BGWilsonSU3_MatMult_NND(res_ptr+str2, 
				U_nu_ptr+str2,UpNu_ptr+str2,UpMu_ptr+str2,ns);
	BGWilsonSU3_MatMult_ND(UpNu_ptr+str2,U_mu_ptr+str2,res_ptr+str2,ns);

#pragma omp for reduction(+:rectF)
	for(int site=0; site<Nvol_; ++site)
	  rectF += UpNu_ptr[CC2*site]+ UpNu_ptr[8+CC2*site]+ UpNu_ptr[16+CC2*site]; ///ReTr

      }
    }
  }

  Communicator::instance()->sync();
  plaqF = Communicator::instance()->reduce_sum(plaqF);
  rectF = Communicator::instance()->reduce_sum(rectF);
 
  int NP = CommonPrms::instance()->NP(); 
  double plaq = plaqF/NC_;
  double rect = rectF/NC_;
  double Hgauge = Params.c_plaq*(Nvol_*NP*NDIM_*(NDIM_-1.0)/2.0 -plaq)
    + Params.c_rect*(Nvol_*NP*NDIM_*(NDIM_-1.0) -rect);
  Hgauge *= Params.beta;
  
  _Message(ACTION_VERB_LEVEL, "    [ActionGaugeRect] H = "<<Hgauge<<"\n");
  _Message(1,"    -- Plaquette = "
	   << plaq/(Nvol_*NP*NDIM_*(NDIM_-1.0)/2.0) << "\n");
  return Hgauge;
}




GaugeField ActionGaugeRect::md_force(){
 using namespace Mapping;
  
  const Staples stpl_;
  SUNmat force_mat;
  GaugeField force;
  GaugeField1D force_pl;
  GaugeField1D force_rect;

  GaugeField1D Cup1;
  GaugeField1D Cup2;
  GaugeField1D Cdn1;
  GaugeField1D Cdn2;

  //check speed
  GaugeField1D res;
  GaugeField1D UpMu, UpNu;
  GaugeField1D U_mu, U_nu;

  GaugePtr* force_rect_ptr = (GaugePtr*)force_rect.data.getaddr(0);
  GaugePtr* force_pl_ptr = (GaugePtr*)force_pl.data.getaddr(0);
  double* UpMu_ptr = UpMu.data.getaddr(0);
  double* UpNu_ptr = UpNu.data.getaddr(0);
  double* res_ptr = res.data.getaddr(0);
  double* U_ptr = u_->data.getaddr(0);
  double* U_nu_ptr =u_->data.getaddr(0);// auxiliary now initialized to a real location
  double* U_mu_ptr =u_->data.getaddr(0);// auxiliary
  double* Cdn1_ptr = Cdn1.data.getaddr(0);
  double* Cdn2_ptr = Cdn2.data.getaddr(0);
  double* Cup2_ptr = Cup2.data.getaddr(0);
  double* Cup1_ptr = Cup1.data.getaddr(0);

  BGQThread_Init();
  
#pragma omp parallel 
  {
    const int nid = omp_get_num_threads();
    const int tid = omp_get_thread_num();
    const int is = tid*Nvol_ / nid;
    const int ie = (tid + 1)*Nvol_ / nid;
    const int ns = ie - is;
    const int CC = NC_*NC_;
    const int str = is*CC;
    const int CC2 = 2*NC_*NC_;
    const int str2 = is*CC2;
    
    for(int mu=0; mu<NDIM_; ++mu){
      BGWilsonSU3_MatZero((double*)force_pl_ptr+str2, ns);
      BGWilsonSU3_MatZero((double*)force_rect_ptr+str2,ns);

      U_mu_ptr = U_ptr+mu*Nvol_*CC2; // assumes contiguity of U_mu     
      for(int nu=0; nu<NDIM_; ++nu){
	if (nu == mu) continue;
	U_nu_ptr = U_ptr+nu*Nvol_*CC2;      

	shiftField(UpNu,U_mu_ptr,nu,Forward());// U_mu(x+nu)
	shiftField(UpMu,U_nu_ptr,mu,Forward());// U_nu(x+mu)

	BGWilsonSU3_MatMult_NND(Cup1_ptr+str2,U_nu_ptr+str2,UpNu_ptr+str2,UpMu_ptr+str2,ns);
	BGWilsonSU3_MatMult_NND(Cup2_ptr+str2,U_mu_ptr+str2,UpMu_ptr+str2,UpNu_ptr+str2,ns);

	BGWilsonSU3_MatMult_DNN(res_ptr+str2, U_nu_ptr+str2,U_mu_ptr+str2,UpMu_ptr+str2,ns);
	BGWilsonSU3_MatMult_DNN(UpMu_ptr+str2,U_mu_ptr+str2,U_nu_ptr+str2,UpNu_ptr+str2,ns);

	shiftField(Cdn1,res_ptr,nu,Backward());
	shiftField(Cdn2,UpMu_ptr,mu,Backward());
	shiftField(UpMu, Cup2_ptr, mu, Forward()); // Cup2(x+mu)
        
 	// plaquette term
	// force_pl += Cup1;
	BGWilsonSU3_MatAdd((double*)force_pl_ptr+str2, 
			   Cup1_ptr+str2, ns); 

	//force_pl += Cdn1;
	BGWilsonSU3_MatAdd((double*)force_pl_ptr+str2, 
			   Cdn1_ptr+str2, ns); 
	
	BGWilsonSU3_MatMult_NND(res_ptr+str2,U_nu_ptr+str2,UpNu_ptr+str2,UpMu_ptr+str2,ns);

	//force_rect += res; 
	BGWilsonSU3_MatAdd((double*)force_rect_ptr+str2, 
			   res_ptr+str2,ns); 

	BGWilsonSU3_MatMult_DNN(Cup2_ptr+str2,U_nu_ptr+str2,U_mu_ptr+str2,UpMu_ptr+str2,ns);
	
	shiftField(UpMu,U_nu_ptr,mu,Forward());    // U_nu(x+mu)
	shiftField(UpNu,Cup1_ptr,nu,Forward());    // Cup1(x+nu)

	BGWilsonSU3_MatMult_NND(res_ptr+str2,U_nu_ptr+str2,UpNu_ptr+str2,UpMu_ptr+str2,ns);

	BGQThread_Barrier(0, nid);
	//force_rect += res; 
	BGWilsonSU3_MatAdd((double*)force_rect_ptr+str2, 
			   res_ptr+str2,ns); 

	BGWilsonSU3_MatMult_DNN(res_ptr+str2,U_nu_ptr+str2,Cdn1_ptr+str2,UpMu_ptr+str2,ns);
	BGWilsonSU3_MatMult_DNN(UpNu_ptr+str2,Cdn2_ptr+str2,U_mu_ptr+str2,UpMu_ptr+str2,ns);

	//res += UpNu;
	BGWilsonSU3_MatAdd(res_ptr+str2, 
			   UpNu_ptr+str2,ns); 
	//res += Cup2;
	BGWilsonSU3_MatAdd(res_ptr+str2, 
			   Cup2_ptr+str2,ns); 

	shiftField(Cup1,res_ptr,nu,Backward());//+=res(x-nu)
	//force_rect += Cup1_ptr; 
	BGWilsonSU3_MatAdd((double*)force_rect_ptr+str2, 
			   Cup1_ptr+str2,ns); 

	shiftField(UpNu,U_mu_ptr,nu,Forward()); // U_mu(x+nu)
	BGWilsonSU3_MatMult_NND(res_ptr+str2,Cdn2_ptr+str2,UpNu_ptr+str2,UpMu_ptr+str2,ns);
	//force_rect += res;
	BGWilsonSU3_MatAdd((double*)force_rect_ptr+str2, 
			   res_ptr+str2,ns); 
      }
      //force_pl   *= Params.c_plaq;
      BGWilsonSU3_MatMultScalar((double*)force_pl_ptr+str2,
				Params.c_plaq,ns);

      //force_rect *= Params.c_rect;
      BGWilsonSU3_MatMultScalar((double*)force_rect_ptr+str2,
				Params.c_rect,ns);
      //force_rect += force_pl; //force_rect = total force (staples term)
      BGWilsonSU3_MatAdd((double*)force_rect_ptr+str2, 
			 (double*)force_pl_ptr+str2,ns); 

      BGWilsonSU3_MatMult_ND(Cdn1_ptr+str2,U_mu_ptr+str2, 
			     force_rect.data.getaddr(0)+str2,ns);


#pragma omp single 
      {    
	SetSlice(force, TracelessAntihermite(Cdn1),mu);
      }
      
    }

  }

  force *= 0.5*Params.beta/NC_;
  _MonitorMsg(ACTION_VERB_LEVEL,Action,force,"ActionGaugeRect");

  return force;
}


