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
  __complex__ double v[9];
}GaugePtr;



double ActionGaugeRect::calc_H(){
  using namespace Mapping;

  const Staples stpl_;
  GaugeField1D Cup1, Cup2;
  GaugeField1D U_nu, U_mu, res;
  GaugeField1D UpNu, UpMu;

  double plaqF = 0.0;
  double rectF = 0.0;

#ifdef IBM_BGQ_WILSON
  //pointers
  double* UpMu_ptr = UpMu.data.getaddr(0);
  double* UpNu_ptr = UpNu.data.getaddr(0);
  double* res_ptr = res.data.getaddr(0);
  double* U_ptr = u_->data.getaddr(0);
  double* U_nu_ptr = U_nu.data.getaddr(0);
  double* U_mu_ptr = U_mu.data.getaddr(0);
  double* Cup1_ptr = Cup1.data.getaddr(0);
  double* Cup2_ptr = Cup2.data.getaddr(0);

  BGQThread_Init();

  // From Matsufuru-san code
#pragma omp parallel
  {
    
    int tid, nid;
    int is, ie, ns;
    nid = omp_get_num_threads();
    tid = omp_get_thread_num();
    is = tid*Nvol_ / nid;
    ie = (tid + 1)*Nvol_ / nid;
    ns = ie - is;
    int jump2 = is*18;
    int jump = is*9;
    
    for(int mu = 0; mu < NDIM_; ++mu){
      BGWilsonLA_MatEquate((__complex__ double*)U_mu_ptr+jump,
			   (__complex__ double*)U_ptr+mu*Nvol_*9+jump, ns); //U_mu links    
      
      
      for(int nu = mu+1; nu < NDIM_; ++nu){
	
	BGWilsonLA_MatEquate((__complex__ double*)U_nu_ptr+jump,
			     (__complex__ double*)U_ptr+jump+nu*Nvol_*9, ns); //U_nu links    
	
	
	shiftField(UpNu,U_mu_ptr,nu,Forward());
	shiftField(UpMu,U_nu_ptr,mu,Forward());
	
	BGWilsonSU3_MatMult_NND(Cup1_ptr+jump2, U_nu_ptr+jump2, UpNu_ptr+jump2, UpMu_ptr+jump2, ns);
	BGWilsonSU3_MatMult_NND(Cup2_ptr+jump2, U_mu_ptr+jump2, UpMu_ptr+jump2, UpNu_ptr+jump2, ns);
	
	// plaquette term
	BGWilsonSU3_MatMult_ND(res_ptr+jump2, U_ptr+mu*Nvol_*18+jump2, Cup1_ptr+jump2, ns);
      
	
#pragma omp for reduction(+:plaqF) //nb NOT "parallel for" just "for"
	for(int s=0; s<Nvol_; s++)
	    plaqF += res_ptr[18*s]+res_ptr[8+18*s]+res_ptr[16+18*s];

	shiftField(UpNu, U_mu_ptr,nu,Forward());      //U_mu(x+nu)
	shiftField(UpMu, Cup2_ptr,mu,Forward());      //Cup2(x+mu) 
	
	BGWilsonSU3_MatMult_NND(res_ptr+jump2, U_nu_ptr+jump2, UpNu_ptr+jump2, UpMu_ptr+jump2, ns);
	BGWilsonSU3_MatMult_ND(UpNu_ptr+jump2, U_ptr+mu*Nvol_*18+jump2, res_ptr+jump2, ns);
	
#pragma omp for reduction(+:rectF)
	for(int site=0; site<Nvol_; ++site)
	    rectF += UpNu_ptr[18*site]+ UpNu_ptr[8+18*site]+ UpNu_ptr[16+18*site];
	
	shiftField(UpNu,Cup1_ptr,nu,Forward());   //Cup1(x+nu)
	shiftField(UpMu,U_nu_ptr,mu,Forward());   //U_nu(x+mu)
	
	BGWilsonSU3_MatMult_NND(res_ptr+jump2, U_nu_ptr+jump2, UpNu_ptr+jump2, UpMu_ptr+jump2, ns);
	BGWilsonSU3_MatMult_ND(UpNu_ptr+jump2, U_ptr+mu*Nvol_*18+jump2, res_ptr+jump2, ns);
	
#pragma omp for reduction(+:rectF)
	  for(int site=0; site<Nvol_; ++site)
	    rectF += UpNu_ptr[18*site]+ UpNu_ptr[8+18*site]+ UpNu_ptr[16+18*site];

      }
    }
  }
#else 
  for(int mu = 0; mu < NDIM_; ++mu){
    U_mu = DirSlice(*u_,mu); 
    
    for(int nu = mu+1; nu < NDIM_; ++nu){
      Cup1 = stpl_.upper(*u_,mu,nu);
      Cup2 = stpl_.upper(*u_,nu,mu);
      
      // plaquette term
      for(int site=0; site<Nvol_; ++site){
	plaqF += ReTr(mat(*u_,site,mu)*mat_dag(Cup1,site) );
      }
      
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
  
  plaqF = Communicator::instance()->reduce_sum(plaqF);
  rectF = Communicator::instance()->reduce_sum(rectF);
 
  int NP = CommonPrms::instance()->NP(); 
  double plaq = plaqF/NC_;
  double rect = rectF/NC_;
  double Hgauge = Params.c_plaq*(Nvol_*NP*NDIM_*(NDIM_-1.0)/2.0 -plaq)
    + Params.c_rect*(Nvol_*NP*NDIM_*(NDIM_-1.0)     -rect);
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
  GaugeField1D U_mu, U_nu, res;
  GaugeField1D UpMu, UpNu;

#ifdef IBM_BGQ_WILSON
  GaugePtr* force_rect_ptr = (GaugePtr*)force_rect.data.getaddr(0);
  double* UpMu_ptr = UpMu.data.getaddr(0);
  double* UpNu_ptr = UpNu.data.getaddr(0);
  double* res_ptr = res.data.getaddr(0);
  double* U_ptr = u_->data.getaddr(0);
  double* U_nu_ptr = U_nu.data.getaddr(0);
  double* U_mu_ptr = U_mu.data.getaddr(0);
  double* Cdn1_ptr = Cdn1.data.getaddr(0);
  double* Cdn2_ptr = Cdn2.data.getaddr(0);
  double* Cup2_ptr = Cup2.data.getaddr(0);
  double* Cup1_ptr = Cup1.data.getaddr(0);

  BGQThread_Init();
  
#pragma omp parallel 
  {
    int tid, nid;
    int is, ie, ns;
    nid = omp_get_num_threads();
    tid = omp_get_thread_num();
    is = tid*Nvol_ / nid;
    ie = (tid + 1)*Nvol_ / nid;
    ns = ie - is;
    int jump = is*18;
    
    for(int mu=0; mu<NDIM_; ++mu){
      BGWilsonLA_MatZero((__complex__ double*)force_pl.data.getaddr(0)+is*9, ns);
      BGWilsonLA_MatZero((__complex__ double*)force_rect_ptr+is*9, ns);
      
      BGWilsonLA_MatEquate((__complex__ double*)U_mu_ptr+is*9,
			   (__complex__ double*)U_ptr+is*9+mu*Nvol_*9, ns); //U_mu links    
   
      for(int nu=0; nu<NDIM_; ++nu){
 	if (nu == mu) continue;
	BGWilsonLA_MatEquate((__complex__ double*)U_nu_ptr+is*9,
			     (__complex__ double*)U_ptr+is*9+nu*Nvol_*9, ns); //U_nu links    
	

	shiftField(UpNu,U_mu_ptr,nu,Forward());// U_mu(x+nu)
	shiftField(UpMu,U_nu_ptr,mu,Forward());// U_nu(x+mu)
	
	BGWilsonSU3_MatMult_NND(Cup1_ptr+jump, U_nu_ptr+jump, UpNu_ptr+jump, UpMu_ptr+jump, ns);
	BGWilsonSU3_MatMult_NND(Cup2_ptr+jump, U_mu_ptr+jump, UpMu_ptr+jump, UpNu_ptr+jump, ns);
	
	BGWilsonSU3_MatMult_DNN(res_ptr+jump, U_nu_ptr+jump, U_mu_ptr+jump, UpMu_ptr+jump, ns);
	BGWilsonSU3_MatMult_DNN(UpMu_ptr+jump, U_mu_ptr+jump, U_nu_ptr+jump, UpNu_ptr+jump, ns);
	
	shiftField(Cdn1,res_ptr,nu,Backward());
	shiftField(Cdn2,UpMu_ptr,mu,Backward());
	shiftField(UpMu, Cup2_ptr, mu, Forward()); // Cup2(x+mu)
        
 	// plaquette term
	// force_pl += Cup1;
	BGWilsonMatLA_Add((__complex__ double*)force_pl.data.getaddr(0)+is*9, 
			  (__complex__ double*)Cup1_ptr+is*9, ns); 
	//force_pl += Cdn1;
	BGWilsonMatLA_Add((__complex__ double*)force_pl.data.getaddr(0)+is*9, 
			  (__complex__ double*)Cdn1_ptr+is*9, ns); 
	
	
	BGWilsonSU3_MatMult_NND(res_ptr+jump, U_nu_ptr+jump, UpNu_ptr+jump, UpMu_ptr+jump, ns);
	
	//force_rect += res;
	BGWilsonMatLA_Add((__complex__ double*)force_rect_ptr+is*9, 
			  (__complex__ double*)res_ptr+is*9, ns); 
	
	
	BGWilsonSU3_MatMult_DNN(Cup2_ptr+jump, U_nu_ptr+jump, U_mu_ptr+jump, UpMu_ptr+jump, ns);   
	
	shiftField(UpMu, U_nu_ptr, mu, Forward());    // U_nu(x+mu)
	shiftField(UpNu, Cup1_ptr, nu, Forward());    // Cup1(x+nu)
	
	BGWilsonSU3_MatMult_NND(res_ptr+jump, U_nu_ptr+jump, UpNu_ptr+jump, UpMu_ptr+jump, ns);

	//force_rect += res; 
	BGWilsonMatLA_Add((__complex__ double*)force_rect_ptr+is*9, 
			  (__complex__ double*)res_ptr+is*9, ns); 

	BGWilsonSU3_MatMult_DNN(res_ptr+jump, U_nu_ptr+jump, Cdn1_ptr+jump, UpMu_ptr+jump, ns);
	BGWilsonSU3_MatMult_DNN(UpNu_ptr+jump, Cdn2_ptr+jump, U_mu_ptr+jump, UpMu_ptr+jump, ns);

	//res += UpNu;
	BGWilsonMatLA_Add((__complex__ double*)res_ptr+is*9, 
			  (__complex__ double*)UpNu_ptr+is*9, ns); 
	//res += Cup2;
	BGWilsonMatLA_Add((__complex__ double*)res_ptr+is*9, 
			  (__complex__ double*)Cup2_ptr+is*9, ns); 

	shiftField(Cup1, res_ptr,nu,Backward());// res(x-nu)
	
	//force_rect += Cup1_ptr;
	BGWilsonMatLA_Add((__complex__ double*)force_rect_ptr+is*9, 
			  (__complex__ double*)Cup1_ptr+is*9, ns); 
	
	shiftField(UpNu, U_mu_ptr, nu, Forward()); // U_mu(x+nu)
	BGWilsonSU3_MatMult_NND(res_ptr+jump, Cdn2_ptr+jump, UpNu_ptr+jump, UpMu_ptr+jump, ns);

	//force_rect += res;
	BGWilsonMatLA_Add((__complex__ double*)force_rect_ptr+is*9, 
			  (__complex__ double*)res_ptr+is*9, ns); 
	
      }
      //force_pl   *= Params.c_plaq;
      BGWilsonLA_MatMultScalar((__complex__ double*)force_pl.data.getaddr(0)+is*9,Params.c_plaq, ns);
      //force_rect *= Params.c_rect;
      BGWilsonLA_MatMultScalar((__complex__ double*)force_rect_ptr+is*9,Params.c_rect, ns);
      
      //force_rect += force_pl; //force_rect = total force (staples term)
      BGWilsonMatLA_Add((__complex__ double*)force_rect_ptr+is*9, 
			(__complex__ double*)force_pl.data.getaddr(0)+is*9, ns); 
       
      BGWilsonSU3_MatMult_ND(Cdn1_ptr+jump , u_->data.getaddr(0)+18*Nvol_*mu+jump, 
			     force_rect.data.getaddr(0)+jump, ns);

#pragma omp single 
      {    
	SetSlice(force, TracelessAntihermite(Cdn1), mu);
      }
      
    }

  }

#else
  for(int mu=0; mu<NDIM_; ++mu){
    force_pl   = 0.0;
    force_rect = 0.0;

    U_mu = DirSlice(*u_, mu);   //U_mu links

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
      U_nu = DirSlice(*u_, nu);   //U_nu links    
      //          U_mu
      //         +-->--+-->--+
      //   U_nu  |           |   term  (Cup2)
      //        (x)    +--<--+      

      UpMu = shiftField(Cup2, mu, Forward()); // Cup2(x+mu)
      UpNu = shiftField(U_mu, nu, Forward()); // U_mu(x+nu)

      res = 0.0;
      for(int site=0; site<Nvol_; ++site){
        res.data[res.format.cslice(0,site)] = 
          (mat(U_nu,site)*mat(UpNu,site)*mat_dag(UpMu,site)).getva();
      }

      force_rect += res;  
  
      //         +-->--+
      //         |     |
      //         +     +   term
      //   U_nu  |     |  U_nu(x+mu) (UpMu)
      //        (x)    v
      UpMu = shiftField(U_nu, mu, Forward());    // U_nu(x+mu)
      UpNu = shiftField(Cup1, nu, Forward());    // Cup1(x+nu)
      res = 0.0;

      for(int site=0; site<Nvol_; ++site){
	res.data[res.format.cslice(0,site)] =
	  (mat(U_nu,site)*mat(UpNu,site)*mat_dag(UpMu,site)).getva();
      }       
     force_rect += res;     
      //           U_mu(x+nu)
      //      +-->--+-->--+
      //      |           |   term
      //      +--<-(x)    v
      UpNu = shiftField(U_mu, nu, Forward()); // U_mu(x+nu)
      res = 0.0; 

      for(int site=0; site<Nvol_; ++site){
	res.data[res.format.cslice(0,site)] = 
	  (mat(Cdn2,site)*mat(UpNu,site)*mat_dag(UpMu,site)).getva();
      }       
      force_rect += res;
      //     (x)    +--<--+
      //      |           |   term
      //      +-->--+-->--+
      UpMu = shiftField(Cup2, mu, Forward());

      res = 0.0;

      for(int site=0; site<Nvol_; ++site){
	res.data[res.format.cslice(0,site)] = 
	  (mat_dag(U_nu,site)*mat(U_mu,site)*mat(UpMu,site)).getva();
      } 
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

      for(int site = 0; site<Nvol_; ++site){
	res.data[res.format.cslice(0,site)] = 
	  (mat_dag(Cdn2,site)*mat(U_mu,site)*mat(UpMu,site)).getva();
      }       
      force_rect += shiftField(res,nu,Backward());//+=res(x-nu)  
    }
    force_pl   *= Params.c_plaq;
    force_rect *= Params.c_rect;
    force_rect += force_pl; //force_rect = total force (staples term)

    for(int site = 0; site<Nvol_; ++site){
      force_mat = (mat(*u_,site,mu)*mat_dag(force_rect,site));
      SetMat(force, anti_hermite_traceless(force_mat), site, mu);
    }
  } 
#endif
  
  force *= 0.5*Params.beta/NC_;
  _MonitorMsg(ACTION_VERB_LEVEL,Action,force,"ActionGaugeRect");

  return force;
}


