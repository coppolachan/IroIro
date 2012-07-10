/*!
 * @file meson_correlator.hpp
 * @brief Class for calculation of generic meson correlator
 */
#include "meson_correlator.hpp"
#include "Communicator/comm_io.hpp"
#include "Main/Geometry/siteIndex.h"

//assuming gamma5 hermiticity
template <typename F>
const std::vector<double> MesonCorrelator::calculate(const prop_t& q1, 
						     const prop_t& q2) {
  //CCIO::cout <<" Contraction to make up pp-correlator\n";
  F fmt;
  int Nt = CommonPrms::instance()->Nt();

  std::vector<double> correl_local(Nt,0.0);
  
  for(int site=0; site<CommonPrms::instance()->Nvol(); ++site){
    int t = SiteIndex::instance()->t(site);//get t from site index
    //Cycle among spinor and color indexes
    for(int s1=0; s1<Nd_; ++s1){
      for(int s2=0; s2<Nd_; ++s2){
	for(int c1=0; c1<Nc_; ++c1){
	  correl_local[t] +=(q1[c1+Nc_*s1][fmt->cslice(s2,site)]
			    *q2[c1+Nc_*s1][fmt->cslice(s2,site)]).sum();
	}
      }
    }
  }
  
  int Lt = CommonPrms::instance()->Lt();
  int ipet = Communicator::instance()->ipe(3);

  std::vector<double> correl_tmp(Lt,0.0);
  for(int t = 0; t< Nt; ++t) correl_tmp[t +ipet*Nt] = correl_local[t];

  std::vector<double> correl(Lt);
  for(int t = 0; t< Lt; ++t) 
    correl[t] = Communicator::instance()->reduce_sum(correl_tmp[t]);
  
  return correl;
}
