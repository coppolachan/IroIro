/*!
 * @file meson_correlator.hpp
 *
 *
 * @brief Class for calculation of generic meson correlator
 *
 */

#include "meson_correlator.hpp"
#include "Communicator/comm_io.hpp"

const std::vector<double> MesonCorrelator::calculate(const prop_t& q1, const prop_t& q2) {
  CCIO::cout <<" Contraction to make up pp-correlator\n";
  
  std::vector<double> correl_local(Nt_,0.0);
  
  for(int site = 0; site < CommonPrms::instance()->Nvol(); ++site){
    int t = SiteIndex::instance()->t(site);//get t from site index
    //Cycle among spinor and color indexes
    for(int spinor1 = 0; spinor1 < Nd_; ++spinor1){
      for(int color1 = 0; color1 < Nc_; ++color1){
	for(int spinor2=0; spinor2 < Nd_; ++spinor2){
	  correl_local[t] +=(q1[color1+Nc_*spinor1][fmt_->cslice(spinor2,site)]
			     *q2[color1+Nc_*spinor1][fmt_->cslice(spinor2,site)]).sum();
	}
      }
    }
  }
  
  int Lt = CommonPrms::instance()->Lt();
  int ipet = Communicator::instance()->ipe(3);

  std::vector<double> correl_tmp(Lt,0.0);
  for(int t = 0; t< Nt_; ++t) correl_tmp[t +ipet*Nt_] = correl_local[t];

  std::vector<double> correl(Lt);
  for(int t = 0; t< Lt; ++t) 
    correl[t] = Communicator::instance()->reduce_sum(correl_tmp[t]);
  
  return correl;



}
