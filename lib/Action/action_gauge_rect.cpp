/*!
  @file action_gauge_rect.cpp
  @brief Definition of the ActionGaugeRect class
*/


#include "action_gauge_rect.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/fieldUtils.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include "lib/Main/Geometry/mapping.hpp"
#include "include/messages_macros.hpp"
#ifdef IBM_BGQ_WILSON
#include <omp.h>
#include "bgqthread.h"
#include "include/bgq_su3algebra.h"
#endif

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

  double plaqF = 0.0;
  double rectF = 0.0;

#ifdef IBM_BGQ_WILSON
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
#else 
  for(int mu=0; mu<NDIM_; ++mu){
    U_mu = DirSlice(*u_,mu); 
    
    for(int nu = mu+1; nu < NDIM_; ++nu){
      Cup1 = stpl_.upper(*u_,mu,nu);
      Cup2 = stpl_.upper(*u_,nu,mu);
      
      // plaquette term
      for(int site=0; site<Nvol_; ++site)
	plaqF += ReTr(mat(*u_,site,mu)*mat_dag(Cup1,site));
      
      U_nu = DirSlice(*u_,nu);
      // rectangular terms
      //       mu,U_mu (UpNu) 
      //         +-->--+-->--+
      // nu,U_nu |           |Cup2_dag(site+mu) (UpMu)
      //     site+     +--<--+                               
      
      UpNu = shiftField(U_mu,nu,Forward());      //U_mu(x+nu)
      UpMu = shiftField(Cup2,mu,Forward());      //Cup2(x+mu)     
      
      res = 0.0;
      for(int site = 0; site<Nvol_; ++site){
	res.data[res.format.cslice(0,site)] = 
	  (mat(U_nu,site)*mat(UpNu,site)*mat_dag(UpMu,site)).getva();
	rectF += ReTr(mat(*u_,site,mu)*mat_dag(res,site));
      }
      
      //        Cup1(site+nu)
      //          +-->--+
      //          |     |
      //          +     +
      // nu,U_nu  |     | U_nu+mu (UpMu)
      //    site  +     +
      
      UpNu = shiftField(Cup1,nu,Forward());   //Cup1(x+nu)
      UpMu = shiftField(U_nu,mu,Forward());   //U_nu(x+mu)
      
      res = 0.0;
      for(int site = 0; site<Nvol_; ++site){
	res.data[res.format.cslice(0,site)] =
	  (mat(U_nu,site)*mat(UpNu,site)*mat_dag(UpMu,site)).getva();
	rectF += ReTr(mat(*u_,site,mu)*mat_dag(res,site));
      }
    }
  }
#endif
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

#ifdef IBM_BGQ_WILSON
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
      BGWilsonLA_MatZero((__complex__ double*)force_pl_ptr+str,ns);
      BGWilsonLA_MatZero((__complex__ double*)force_rect_ptr+str,ns);

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
	BGWilsonMatLA_Add((__complex__ double*)force_pl_ptr+str, 
			  (__complex__ double*)Cup1_ptr+str, ns); 

	//force_pl += Cdn1;
	BGWilsonMatLA_Add((__complex__ double*)force_pl_ptr+str, 
			  (__complex__ double*)Cdn1_ptr+str, ns); 
	
	BGWilsonSU3_MatMult_NND(res_ptr+str2,U_nu_ptr+str2,UpNu_ptr+str2,UpMu_ptr+str2,ns);

	//force_rect += res; 
	BGWilsonMatLA_Add((__complex__ double*)force_rect_ptr+str, 
			  (__complex__ double*)res_ptr+str,ns); 

	BGWilsonSU3_MatMult_DNN(Cup2_ptr+str2,U_nu_ptr+str2,U_mu_ptr+str2,UpMu_ptr+str2,ns);
	
	shiftField(UpMu,U_nu_ptr,mu,Forward());    // U_nu(x+mu)
	shiftField(UpNu,Cup1_ptr,nu,Forward());    // Cup1(x+nu)

	BGWilsonSU3_MatMult_NND(res_ptr+str2,U_nu_ptr+str2,UpNu_ptr+str2,UpMu_ptr+str2,ns);

	BGQThread_Barrier(0, nid);
	//force_rect += res; 
	BGWilsonMatLA_Add((__complex__ double*)force_rect_ptr+str, 
			  (__complex__ double*)res_ptr+str,ns); 

	BGWilsonSU3_MatMult_DNN(res_ptr+str2,U_nu_ptr+str2,Cdn1_ptr+str2,UpMu_ptr+str2,ns);
	BGWilsonSU3_MatMult_DNN(UpNu_ptr+str2,Cdn2_ptr+str2,U_mu_ptr+str2,UpMu_ptr+str2,ns);

	//res += UpNu;
	BGWilsonMatLA_Add((__complex__ double*)res_ptr+str, 
			  (__complex__ double*)UpNu_ptr+str,ns); 
	//res += Cup2;
	BGWilsonMatLA_Add((__complex__ double*)res_ptr+str, 
			  (__complex__ double*)Cup2_ptr+str,ns); 

	shiftField(Cup1,res_ptr,nu,Backward());//+=res(x-nu)
	//force_rect += Cup1_ptr; 
	BGWilsonMatLA_Add((__complex__ double*)force_rect_ptr+str, 
			  (__complex__ double*)Cup1_ptr+str,ns); 

	shiftField(UpNu,U_mu_ptr,nu,Forward()); // U_mu(x+nu)
	BGWilsonSU3_MatMult_NND(res_ptr+str2,Cdn2_ptr+str2,UpNu_ptr+str2,UpMu_ptr+str2,ns);
	//force_rect += res;
	BGWilsonMatLA_Add((__complex__ double*)force_rect_ptr+str, 
			  (__complex__ double*)res_ptr+str,ns); 
      }
      //force_pl   *= Params.c_plaq;
      BGWilsonLA_MatMultScalar((__complex__ double*)force_pl_ptr+str,
			       Params.c_plaq,ns);

      //force_rect *= Params.c_rect;
      BGWilsonLA_MatMultScalar((__complex__ double*)force_rect_ptr+str,
			       Params.c_rect,ns);
      //force_rect += force_pl; //force_rect = total force (staples term)
      BGWilsonMatLA_Add((__complex__ double*)force_rect_ptr+str, 
			(__complex__ double*)force_pl_ptr+str,ns); 

      BGWilsonSU3_MatMult_ND(Cdn1_ptr+str2,U_mu_ptr+str2, 
			     force_rect.data.getaddr(0)+str2,ns);


#pragma omp single 
      {    
	SetSlice(force, TracelessAntihermite(Cdn1),mu);
      }
      
    }

  }

#else
  for(int mu=0; mu<NDIM_; ++mu){
    force_pl   = 0.0;
    force_rect = 0.0;

    U_mu = DirSlice(*u_,mu);   //U_mu links

    for(int nu=0; nu<NDIM_; ++nu){
      if (nu == mu) continue;
      
      Cup1 = stpl_.upper(*u_,mu,nu);
      Cup2 = stpl_.upper(*u_,nu,mu);
      Cdn1 = stpl_.lower(*u_,mu,nu);
      Cdn2 = stpl_.lower(*u_,nu,mu);

      // plaquette term
      force_pl += Cup1;
      force_pl += Cdn1;

      // rectangular terms
      // ^nu
      // |  
      // +-->mu
      //
      // (x) is the site position
      U_nu = DirSlice(*u_,nu);   //U_nu links    
      //          U_mu
      //         +-->--+-->--+
      //   U_nu  |           |   term  (Cup2)
      //        (x)    +--<--+      

      UpMu = shiftField(Cup2,mu,Forward()); // Cup2(x+mu)
      UpNu = shiftField(U_mu,nu,Forward()); // U_mu(x+nu)

      res = 0.0;
      for(int site=0; site<Nvol_; ++site)
        res.data[res.format.cslice(0,site)] = 
          (mat(U_nu,site)*mat(UpNu,site)*mat_dag(UpMu,site)).getva();

      force_rect += res;  
  
      //         +-->--+
      //         |     |
      //         +     +   term
      //   U_nu  |     |  U_nu(x+mu) (UpMu)
      //        (x)    v
      UpMu = shiftField(U_nu, mu, Forward());    // U_nu(x+mu)
      UpNu = shiftField(Cup1, nu, Forward());    // Cup1(x+nu)
      res = 0.0;

      for(int site=0; site<Nvol_; ++site)
	res.data[res.format.cslice(0,site)] =
	  (mat(U_nu,site)*mat(UpNu,site)*mat_dag(UpMu,site)).getva();
      
     force_rect += res;     
      //           U_mu(x+nu)
      //      +-->--+-->--+
      //      |           |   term
      //      +--<-(x)    v
      UpNu = shiftField(U_mu,nu,Forward()); // U_mu(x+nu)

      res = 0.0; 
      for(int site=0; site<Nvol_; ++site)
	res.data[res.format.cslice(0,site)] = 
	  (mat(Cdn2,site)*mat(UpNu,site)*mat_dag(UpMu,site)).getva();

      force_rect += res;
      //     (x)    +--<--+
      //      |           |   term
      //      +-->--+-->--+
      UpMu = shiftField(Cup2,mu,Forward());

      res = 0.0;
      for(int site=0; site<Nvol_; ++site)
	res.data[res.format.cslice(0,site)] = 
	  (mat_dag(U_nu,site)*mat(U_mu,site)*mat(UpMu,site)).getva();

      force_rect += shiftField(res,nu,Backward()); //+=res(x-nu)
      //     (x)    ^
      //      |     |
      //      +     +   term
      //      |     |
      //      +-->--+
      UpMu = shiftField(U_nu, mu, Forward()); //U_nu(x+mu)

      res = 0.0;
      for(int site=0; site<Nvol_; ++site){
	res.data[res.format.cslice(0,site)] = 
	  (mat_dag(U_nu,site)*mat(Cdn1,site)*mat(UpMu,site)).getva();
      } 
      force_rect += shiftField(res,nu,Backward());//+=res(x-nu)    
      //      +--<-(x)    ^
      //      |           |   term
      //      +--<--+--<--+

      res = 0.0;
      for(int site=0; site<Nvol_; ++site)
	res.data[res.format.cslice(0,site)] = 
	  (mat_dag(Cdn2,site)*mat(U_mu,site)*mat(UpMu,site)).getva();

      force_rect += shiftField(res,nu,Backward());//+=res(x-nu)  
    }
    force_pl   *= Params.c_plaq;
    force_rect *= Params.c_rect;
    force_rect += force_pl; //force_rect = total force (staples term)

    for(int site=0; site<Nvol_; ++site){
      force_mat = (mat(*u_,site,mu)*mat_dag(force_rect,site));
      SetMat(force, anti_hermite_traceless(force_mat), site, mu);
    }
  } 
#endif  
  force *= 0.5*Params.beta/NC_;
  _MonitorMsg(ACTION_VERB_LEVEL,Action,force,"ActionGaugeRect");

  return force;
}


