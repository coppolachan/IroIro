#include "meson_correlator.hpp"
#include "include/field.h"
#include "Communicator/comm_io.hpp"
#include "Main/Geometry/siteIndex.hpp"
const std::vector<double> 
MesonCorrelator::calculate(const prop_t& q1,const prop_t& q2){

  Format::Format_F fmt(CommonPrms::instance()->Nvol());
  int Nt = CommonPrms::instance()->Nt();

  std::vector<double> correl_local(Nt,0.0);
  int s1, s2, s3, s4; //spinor indexes
  double gamma_factor, temp_corr;

  for(int site=0; site<CommonPrms::instance()->Nvol(); ++site){
    int t = SiteIndex::instance()->c_t(site);//get t from site index
    //Cycle among spinor and color indexes
    for(s4=0; s4<Nd_; ++s4){
      // s4 -> s1
      s1 = (*G2_)(s4).spinor_index; 
      for(s2=0; s2<Nd_; ++s2){
	// s2 -> s3
	s3 = (*G1_)(s2).spinor_index; 
	gamma_factor 
	  = (*G2_)(s4).complex_factor[0]*(*G1_)(s2).complex_factor[0]
	   -(*G2_)(s4).complex_factor[1]*(*G1_)(s2).complex_factor[1]; 
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
  int ipet = Communicator::instance()->ipe(3);

  std::vector<double> correl_tmp(Lt,0.0);
  for(int t = 0; t< Nt; ++t) correl_tmp[t +ipet*Nt] = correl_local[t];

  std::vector<double> correl(Lt);
  for(int t = 0; t< Lt; ++t) 
    correl[t] = Communicator::instance()->reduce_sum(correl_tmp[t]);
  
  return correl;
}
