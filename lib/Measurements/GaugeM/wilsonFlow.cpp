/*!@file wilsonFlow.cpp
 * @brief implementation of the WilsonFlow class
 */
#include "wilsonFlow.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/fieldUtils.hpp"
#include "Communicator/fields_io.hpp"
#ifdef IBM_BGQ_WILSON
#include <omp.h>
#include "bgqthread.h"
#endif
using namespace std;

void WilsonFlow::update_U(const GaugeField& Z)const{
  using namespace SUNmatUtils;
  using namespace FieldUtils;

  double eps = -2.0*estep_;
  int Nvol = CommonPrms::instance()->Nvol();

  for(int m=0; m<NDIM_; ++m){
    for(int site=0; site<Nvol; ++site){
      SUNmat au = exponential(mat(Z,site,m)*eps,Nexp_);
      au.reunit();
      au *= mat(U_,site,m);
      U_.data.set(U_.format.islice(site,m),au.getva());
    }
  }
}

void WilsonFlow::evolve_step()const{
                                   // W0 = U_
  GaugeField Z = Sg_->md_force();  // Z0
  Z *= 0.25;
  update_U(Z);                     // U_= W1 = exp(ep*Z0/4)*W0 

  Z *= -17.0/8.0;
  Z += Sg_->md_force();            // -17/32*Z0 +Z1 
  Z *= 8.0/9.0;                    // -17/36*Z0 +8/9*Z1
  update_U(Z);                     // U_= W2 = exp(ep*(-17/36*Z0 +8/9*Z1))*W1
  
  Z *= -4.0/3.0;
  Z += Sg_->md_force();            // 4/3*(17/36*Z0 -8/9*Z1) +Z2       
  Z *= 0.75;                       // 17/36*Z0 -8/9*Z1 +3/4*Z2       
  update_U(Z);              // V(t+e) = exp(ep*(17/36*Z0 -8/9*Z1 +3/4*Z2))*W2
}
/*
void WilsonFlow::evolve_step()const{
                                    // W0 = U_
  GaugeField Z0 = Sg_->md_force();  // Z0
  Z0 *= 0.5;
  update_U(Z0);                     // U_= W1 = exp(ep*Z0/2)*W0 
  Z0 *= 2.0;

  GaugeField Z1 = Sg_->md_force();  // Z1
  GaugeField Zt = Z1;
  Zt -= Z0;
  Zt *= 1.0-sqrt(2.0)*0.5;          // (1-1/sqrt(2))*(Z1-Z0)  
  update_U(Zt);                     // U_= W2 = exp(ep*Z0/2)*W1 

  GaugeField Z2 = Sg_->md_force();  // Z2
  Z2 *= 1.0+sqrt(2.0)*0.5;          // (1+1/sqrt(2))*Z2
  Zt = Z2;
  Zt -= Z1;
  Z0 *= -0.5*(sqrt(2.0)-1.0);
  Zt += Z0;                        // (1+1/sqrt(2))*Z2 -Z1 +(1/2-1/sqrt(2))*Z0
  Z0 *= -2.0*(sqrt(2.0)+1.0);       
  update_U(Zt);// U_= W3= exp(ep*((1+1/sqrt(2))*Z2 -Z1 +(1/2-1/sqrt(2))*Z0))*W1 

  Zt = Sg_->md_force();     // Z3       
  Z2 *= -4.0;             
  Z1 *= 2.0*(1.0 +sqrt(2.0)); 
  Zt += Z2;                 // Z3 -4*c2*Z2 
  Zt += Z1;                 // Z3 -4*c2*Z2 +(6-4*c1)*Z1
  Zt += Z0;                 // Z3 -4*c2*Z2 +(6-4*c1)*Z1 +Z0
  Zt /= 6.0;
  update_U(Zt);          
}
*/

double WilsonFlow::Edens_plaq(int t)const{
  double td = tau(t);
  return 2.0*td*td*Sg_->calc_H()/double(CommonPrms::instance()->Lvol());
}

double WilsonFlow::Edens_clover(int t)const{
  int Ndim = CommonPrms::instance()->Ndim();
  int Nvol = CommonPrms::instance()->Nvol();

  Staples stpl;
  double Sg = 0.0;      

#ifdef IBM_BGQ_WILSON
  GaugeField1D Fmn;
  double* Fmn_ptr = Fmn.data.getaddr(0);

  BGQThread_Init();

  const int CC2 = 2*NC_*NC_;
  for(int mu=0; mu<Ndim; ++mu){
    for(int nu=mu+1; nu<Ndim; ++nu){
      Fmn = stpl.fieldStrength(U_,mu,nu);

#pragma omp for reduction(+:Sg)	
      for(int site=0; site<Nvol; ++site){
	for (int icc2=0; icc2< CC2; ++icc2){
	  Sg += Fmn_ptr[icc2+CC2*site]*Fmn_ptr[icc2+CC2*site];
	}
      }
    }
  }
#else
  using namespace SUNmatUtils;
  using namespace FieldUtils;

  for(int mu=0; mu<Ndim; ++mu){
    for(int nu=mu+1; nu<Ndim; ++nu){
      GaugeField1D Fmn = stpl.fieldStrength(U_,mu,nu);
      for(int site=0; site<Nvol; ++site){
	valarray<double> u = Fmn.data[Fmn.format.islice(site)];
	Sg += (u*u).sum();   // this is equivarent to Sg += Tr(u*u);
      }
    }
  }
#endif
  Sg = Communicator::instance()->reduce_sum(Sg);

  double td = tau(t);
  return td*td*Sg/double(CommonPrms::instance()->Lvol());
}

void WilsonFlow::save_config(const std::string& fname)const{
  if(saveConf_){
    CCIO::cout << "Saving configuration on disk in binary format\n";
    CCIO::SaveOnDisk<Format::Format_G>(U_.data,fname.c_str());
  }
}
