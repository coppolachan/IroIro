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
      SUNmat au = exponential(mat(Z,site,m)*estep_,Nexp_);
      au *= mat(U_,site,m);
      U_.data.set(U_.format.islice(site,m),au.reunit().getva());
    }
  }
}

void WilsonFlow::gradFlow()const{
  //                                  W0 = U_
  GaugeField Z = Sg_->md_force();  // Z0
  Z *= 0.25;
  update_U(Z);                     // W1 = exp(Z0/4)*W0 

  Z *= -17.0/32.0;                    
  Z += Sg_->md_force();            // +Z1
  Z *= 8.0/9.0;
  update_U(Z);                     // W2 = exp(8/9*Z1 -17/36*Z0)*W1

  Z /= -0.75;                 
  Z += Sg_->md_force();            // +Z2
  Z *= 0.75;               
  update_U(Z);                // V(t+e) = exp(3/4*Z2 -8/9*Z1 +17/36*Z0)*W2
}
  
std::vector<double> WilsonFlow::evolve()const{

  std::vector<double> Eplaq(Nstep_);
  for(int t=0; t< Nstep_; ++t){
    gradFlow();
     Eplaq[t] = 2.0/double(Lvol_)*Sg_->calc_H();
  }
  return Eplaq;
}
