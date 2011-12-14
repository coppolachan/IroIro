/*!
 * @file meson_correlator.hpp
 *
 *
 * @brief Class for calculation of generic meson correlator
 *
 */

#ifndef MESON_CORRELATOR_HPP_
#define MESON_CORRELATOR_HPP_

#include <vector>
#include "include/field.h"

#include "meson_correlator.hpp"
#include "Communicator/comm_io.hpp"
#include "Main/Geometry/siteIndex.h"


typedef std::vector<Field> prop_t;
 
class MesonCorrelator {
  int Nc_;
  int Nd_;  
public:
  MesonCorrelator():Nc_(CommonPrms::instance()->Nc()),
		    Nd_(CommonPrms::instance()->Nd()){}
  template <typename F> const std::vector<double> calculate(const prop_t&, const prop_t&);

};


//assuming gamma5 hermiticity
template <typename F>
const std::vector<double> MesonCorrelator::calculate(const prop_t& q1, const prop_t& q2) {
  CCIO::cout <<" Contraction to make up pp-correlator\n";
  F fmt(CommonPrms::instance()->Nvol());
  int Nt = CommonPrms::instance()->Nt();

  std::vector<double> correl_local(Nt,0.0);
  
  for(int site = 0; site < CommonPrms::instance()->Nvol(); ++site){
    int t = SiteIndex::instance()->t(site);//get t from site index
    //Cycle among spinor and color indexes
    for(int spinor1 = 0; spinor1 < Nd_; ++spinor1){
      for(int color1 = 0; color1 < Nc_; ++color1){
	for(int spinor2=0; spinor2 < Nd_; ++spinor2){
	  correl_local[t] +=(q1[color1+Nc_*spinor1][fmt.cslice(spinor2,site)]
			     *q2[color1+Nc_*spinor1][fmt.cslice(spinor2,site)]).sum();
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

#endif //MESON_CORRELATOR_HPP_
