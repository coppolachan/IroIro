/*!
 * @file meson_correlator.hpp
 * @brief Class for calculation of generic meson correlator
 * It uses the GammaMatrices namespace
 */
#ifndef MESON_CORRELATOR_HPP_
#define MESON_CORRELATOR_HPP_

#include <vector>
#include <memory>

#include "include/field.h"
#include "meson_correlator.hpp"
#include "Communicator/comm_io.hpp"
#include "Main/Geometry/siteIndex.hpp"
#include "Tools/gammaMatrices.hpp"

typedef std::vector<Field> prop_t;
 
enum MesonType {Pion, Scalar, Vector1, Vector2, Vector3};

class MesonCorrelator {
  const int Nc_;
  const int Nd_; 
  const int Nvol_; 
  std::auto_ptr<GammaMatrices::Gamma> G1_;
  std::auto_ptr<GammaMatrices::Gamma> G2_;
public:
  MesonCorrelator(GammaMatrices::Gamma& G1, GammaMatrices::Gamma& G2)
 :Nc_(CommonPrms::instance()->Nc()),
  Nd_(CommonPrms::instance()->Nd()),
  Nvol_(CommonPrms::instance()->Nvol()),
  G1_(&G1), G2_(&G2){}

  MesonCorrelator(MesonType Type = Pion):Nc_(CommonPrms::instance()->Nc()),
					 Nd_(CommonPrms::instance()->Nd()),
					 Nvol_(CommonPrms::instance()->Nvol()){
    switch (Type) {
    case Pion:
      G1_.reset(new GammaMatrices::Unit);
      G2_.reset(new GammaMatrices::Unit);
      break;
    case Scalar:
      G1_.reset(new GammaMatrices::Gamma5);
      G2_.reset(new GammaMatrices::Gamma5);
      break;
    case Vector1:
      G1_.reset(new GammaMatrices::Gamma1_5);
      G2_.reset(new GammaMatrices::Gamma1_5);
      break;
    case Vector2:
      G1_.reset(new GammaMatrices::Gamma2_5);
      G2_.reset(new GammaMatrices::Gamma2_5);
      break;
    case Vector3:
      G1_.reset(new GammaMatrices::Gamma3_5);
      G2_.reset(new GammaMatrices::Gamma3_5);
      break;
    }
  }
  template <typename F> const std::vector<double> calculate(const prop_t&, 
							    const prop_t&); 
};

/*!
 * @brief Calculates the meson correlators with general gamma matrices
 *
 * Using the following structure
 * \f[ \sum_{s_2,s_4,c_1,c_2} G^{s_1,s_2}_{c_1,c_2}(0,x) \Gamma^{s_2,s_3}\cdot G^{{*}{s_4,s_3}}_{c_1,c_2}(0,x) \Gamma^{s_4,s_1}  \f]
 *
 * where $s_1$ is a function of $s_4$ as well as $s_3(s_2)$.
 */
template <typename F>
const std::vector<double> MesonCorrelator::calculate(const prop_t& q1, 
						     const prop_t& q2){
  //CCIO::cout <<"Contraction to make up meson correlator\n";

  F fmt(CommonPrms::instance()->Nvol());
  int Nt = CommonPrms::instance()->Nt();

  std::vector<double> correl_local(Nt,0.0);
  int s1, s2, s3, s4; //spinor indexes
  double gamma_factor, temp_corr;

  for(int site=0; site<CommonPrms::instance()->Nvol(); ++site){
    int t = SiteIndex::instance()->c_t(site);//get t from site index
    //Cycle among spinor and color indexes
    for(s4=0; s4<Nd_; ++s4){
      // s4 -> s1
      s1 = (*G2_)(s4).spn;
      for(s2=0; s2<Nd_; ++s2){
	// s2 -> s3
	s3 = (*G1_)(s2).spn;
	gamma_factor = (*G2_)(s4).facr*(*G1_)(s2).facr
               	      -(*G2_)(s4).faci*(*G1_)(s2).faci; 
	//is always a real number
	for(int c1=0; c1<Nc_; ++c1){//Contracts colors
	  // following contraction is Re[q1*(q2~)] 
	  //  (imaginary part should be exactly zero)
	  temp_corr = (q1[c1+Nc_*s1][fmt.cslice(s2,site)]
		       *q2[c1+Nc_*s4][fmt.cslice(s3,site)]).sum();
	  correl_local[t] += gamma_factor*temp_corr;//no
	}
      }
    }
  }
  
  int Lt = CommonPrms::instance()->Lt();

  std::vector<double> correl_tmp(Lt,0.0);
  for(int t = 0; t< Nt; ++t) 
    correl_tmp[SiteIndex::instance()->global_t(t)] = correl_local[t];

  std::vector<double> correl(Lt);
  for(int t = 0; t< Lt; ++t) 
    correl[t] = Communicator::instance()->reduce_sum(correl_tmp[t]);
  
  return correl;
}


#endif //MESON_CORRELATOR_HPP_
