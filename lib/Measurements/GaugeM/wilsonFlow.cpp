/*!@file wilsonFlow.cpp
 * @brief implementation of the WilsonFlow class
 */

#include "Tools/sunMatUtils.hpp"
#include "Tools/fieldUtils.hpp"
#include "wilsonFlow.hpp"

void WilsonFlow::update_U(const GaugeField& Z)const{
  using namespace SUNmatUtils;
  using namespace FieldUtils;

  for(int m=0; m<NDIM_; ++m){
    for(int site=0; site<Nvol_; ++site){
      SUNmat au = exponential(mat(Z,site,m)*estep_*g0q_,Nexp_);
      au *= mat(U_,site,m);
      U_.data.set(U_.format.islice(site,m),au.reunit().getva());
    }
  }
}

void WilsonFlow::gradFlow()const{
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
void WilsonFlow::gradFlow()const{
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

std::vector<double> WilsonFlow::evolve()const{

  std::vector<double> ttE(Nstep_);
  for(int t=0; t< Nstep_; ++t){
    gradFlow();
    double tau = estep_*(t+1);
    ttE[t] = tau*tau*fabs(g0q_)*Sg_->calc_H()/double(Lvol_);
  }
  return ttE;
}
