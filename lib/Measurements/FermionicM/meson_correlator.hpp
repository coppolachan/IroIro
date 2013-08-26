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
#include "Communicator/comm_io.hpp"
#include "Geometry/siteIndex.hpp"
#include "Tools/gammaMatrices.hpp"
#include "Tools/sunVec.hpp"

typedef std::vector<Field> prop_t;
 
enum MesonType {Pion, Scalar,
		Vector1, Vector2, Vector3, Vector4,
		AVector1, AVector2, AVector3, AVector4,
		Tensor12, Tensor13, Tensor14, Tensor23, Tensor24, Tensor34};

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
    case Vector4:
      G1_.reset(new GammaMatrices::Gamma4_5);
      G2_.reset(new GammaMatrices::Gamma4_5);
      break;
    case AVector1:
      G1_.reset(new GammaMatrices::Gamma1);
      G2_.reset(new GammaMatrices::Gamma1);
      break;
    case AVector2:
      G1_.reset(new GammaMatrices::Gamma2);
      G2_.reset(new GammaMatrices::Gamma2);
      break;
    case AVector3:
      G1_.reset(new GammaMatrices::Gamma3);
      G2_.reset(new GammaMatrices::Gamma3);
      break;
    case AVector4:
      G1_.reset(new GammaMatrices::Gamma4);
      G2_.reset(new GammaMatrices::Gamma4);
      break;
    case Tensor12:
      G1_.reset(new GammaMatrices::Sigma12_5);
      G2_.reset(new GammaMatrices::Sigma12_5);
      break;
    case Tensor13:
      G1_.reset(new GammaMatrices::Sigma13_5);
      G2_.reset(new GammaMatrices::Sigma13_5);
      break;
    case Tensor14:
      G1_.reset(new GammaMatrices::Sigma14_5);
      G2_.reset(new GammaMatrices::Sigma14_5);
      break;
    case Tensor23:
      G1_.reset(new GammaMatrices::Sigma23_5);
      G2_.reset(new GammaMatrices::Sigma23_5);
      break;
    case Tensor24:
      G1_.reset(new GammaMatrices::Sigma24_5);
      G2_.reset(new GammaMatrices::Sigma24_5);
      break;
    case Tensor34:
      G1_.reset(new GammaMatrices::Sigma34_5);
      G2_.reset(new GammaMatrices::Sigma34_5);
      break;
    }
  }
  template <typename F> 
  const std::vector<double> calculate(const prop_t&,const prop_t&);

  template <typename F> 
  const std::vector<double> calculate_mom(const prop_t&,const prop_t&,
					  double a,double b,double c); 
  template <typename F> 
  const std::vector<double> calculate_PSmom(const prop_t&,const prop_t&,
					    double a,double b,double c); 
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

  for(int site=0; site<CommonPrms::instance()->Nvol(); ++site){
    int t = SiteIndex::instance()->c_t(site);//get t from site index
    /// loop over spinor and color indexes
    for(int s4=0; s4<Nd_; ++s4){
      /// s4 -> s1
      int s1 = (*G2_)(s4).spn;
      for(int s2=0; s2<Nd_; ++s2){
        /// s2 -> s3
	int s3 = (*G1_)(s2).spn;
        double gamma_factor
	  = (*G2_)(s4).facr*(*G1_)(s2).facr -(*G2_)(s4).faci*(*G1_)(s2).faci;
        /// is always a real number

	for(int c1=0; c1<Nc_; ++c1){//Contracts colors                        
	  // following contraction is Re[q1*(q2~)]
	  //  (imaginary part should be exactly zero)
          double temp_corr = (q1[c1+Nc_*s1][fmt.cslice(s2,site)]
			      *q2[c1+Nc_*s4][fmt.cslice(s3,site)]).sum();
	  correl_local[t] += gamma_factor*temp_corr;
	}
      }
    }
  }
  int Lt = CommonPrms::instance()->Lt();

  std::vector<double> correl_tmp(Lt,0.0);
  for(int t=0; t<Nt; ++t)
    correl_tmp[SiteIndex::instance()->global_t(t)] = correl_local[t];

  std::vector<double> correl(Lt);
  for(int t=0; t<Lt; ++t)
    correl[t] = Communicator::instance()->reduce_sum(correl_tmp[t]);

  return correl;
}

template <typename F>
const std::vector<double> 
MesonCorrelator::calculate_mom(const prop_t& q1,const prop_t& q2,
			       double px,double py,double pz){

  //CCIO::cout <<"Contraction to make up meson correlator\n";

  F fmt(CommonPrms::instance()->Nvol());
  int Nt = CommonPrms::instance()->Nt();
  std::vector<double> correl_local_re(Nt,0.0),correl_local_im(Nt,0.0);

  for(int site=0; site<CommonPrms::instance()->Nvol(); ++site){
    int t = SiteIndex::instance()->c_t(site);/// get t from site index
    
    /// get global site and space coordinate
    int gsite = SiteIndex::instance()->get_gsite(site);
    int x = SiteIndex::instance()->g_x(gsite);
    int y = SiteIndex::instance()->g_y(gsite);
    int z = SiteIndex::instance()->g_z(gsite);

    /// phase factor
    double sld = px*x + py*y + pz*z;
    double re_ph = cos(sld);
    double im_ph = sin(sld);
	  
    /// loop over spinor and color indices
    for(int s4=0; s4<Nd_; ++s4){
      /// s4 -> s1
      int s1 = (*G2_)(s4).spn;
      for(int s2=0; s2<Nd_; ++s2){
	/// s2 -> s3
	int s3 = (*G1_)(s2).spn;
	double gamma_factor 
	  = (*G2_)(s4).facr*(*G1_)(s2).facr -(*G2_)(s4).faci*(*G1_)(s2).faci; 
	/// is always a real number

	for(int c1=0; c1<Nc_; ++c1){// color contraction
	  /// real part
	  double temp_corr_re = (q1[c1+Nc_*s1][fmt.cslice(s2,site)]
				 *q2[c1+Nc_*s4][fmt.cslice(s3,site)]).sum();
	  /// imaginary part
	  double temp_corr_im = SUNvec(q1[c1+Nc_*s1][fmt.cslice(s2,site)]).im_prod(SUNvec(q2[c1+Nc_*s4][fmt.cslice(s3,site)]));
	  
	  /// multiply the phase
	  correl_local_re[t] += gamma_factor*(temp_corr_re*re_ph - temp_corr_im*im_ph);
	  correl_local_im[t] += gamma_factor*(temp_corr_re*im_ph + temp_corr_im*re_ph);
	}
      }
    }
  }
  int Lt = CommonPrms::instance()->Lt();

  std::vector<double> correl_tmp(2*Lt,0.0);
  for(int t=0; t<Nt; ++t) {
    correl_tmp[2*SiteIndex::instance()->global_t(t)] = correl_local_re[t];
    correl_tmp[2*SiteIndex::instance()->global_t(t)+1] = correl_local_im[t];
  }

  std::vector<double> correl(2*Lt);
  for(int t=0; t<2*Lt; ++t) 
    correl[t] = Communicator::instance()->reduce_sum(correl_tmp[t]);
  
  return correl;
}

template <typename F>
const std::vector<double> 
MesonCorrelator::calculate_PSmom(const prop_t& q1,const prop_t& q2,
			       double px,double py,double pz){

  //CCIO::cout <<"Contraction to make up meson correlator\n";

  F fmt(CommonPrms::instance()->Nvol());
  int Nt = CommonPrms::instance()->Nt();
  std::vector<double> correl_local_re(Nt,0.0),correl_local_im(Nt,0.0);

  for(int site=0; site<CommonPrms::instance()->Nvol(); ++site){
    int t = SiteIndex::instance()->c_t(site);//get t from site index
	  
    /// get global site and space coordinate
    int gsite = SiteIndex::instance()->get_gsite(site);
    int x = SiteIndex::instance()->g_x(gsite);
    int y = SiteIndex::instance()->g_y(gsite);
    int z = SiteIndex::instance()->g_z(gsite);
    
    /// phase factor
    double re_ph = cos(px*x)*cos(py*y)*cos(pz*z);
    
    //Cycle among spinor and color indexes
    for(int s4=0; s4<Nd_; ++s4){
      /// s4 -> s1
      int s1 = (*G2_)(s4).spn;
      for(int s2=0; s2<Nd_; ++s2){
	/// s2 -> s3
	int s3 = (*G1_)(s2).spn;
	double gamma_factor = (*G2_)(s4).facr*(*G1_)(s2).facr
	  -(*G2_)(s4).faci*(*G1_)(s2).faci; 
	/// is always a real number

	for(int c1=0; c1<Nc_; ++c1){//Contracts colors
	  
	  ///real part
	  double temp_corr_re = (q1[c1+Nc_*s1][fmt.cslice(s2,site)]
				 *q2[c1+Nc_*s4][fmt.cslice(s3,site)]).sum();
	  /// multiply the plane wave
	  correl_local_re[t] += gamma_factor*temp_corr_re*re_ph;
	}
      }
    }
  }
  int Lt = CommonPrms::instance()->Lt();

  std::vector<double> correl_tmp(Lt,0.0);
  for(int t=0; t<Nt; ++t)
    correl_tmp[SiteIndex::instance()->global_t(t)] = correl_local_re[t];

  std::vector<double> correl(Lt);
  for(int t=0; t<Lt; ++t) 
    correl[t] = Communicator::instance()->reduce_sum(correl_tmp[t]);
  
  return correl;
}

#endif //MESON_CORRELATOR_HPP_
