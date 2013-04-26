/*! 
  @file mdExec.cpp
  @brief utilities for MD including funcs to generate initial HMC momentum
 */
#include "Geometry/siteIndex.hpp"
#include "mdExec.hpp"
#include "include/format_A.h"
#include "include/format_G.h"
#include "include/field.h"
#include "Tools/randNum_MP.h"
#include "include/timings.hpp"

using namespace std;


void MDexec::md_mom(GaugeField& P,const RandNum& rand){
  if(CommonPrms::instance()->Nc() == 3){
    MDutils::md_mom_su3(P,rand);
  }else{
    MDutils::md_mom_suN(Representation, P,rand);
  }
}

namespace MDutils{
  void md_mom_su3(GaugeField& P,const RandNum& rand){

    static const double sq3i = 1.0/sqrt(3.0);
    int Nvol = CommonPrms::instance()->Nvol();
    int Ndim = CommonPrms::instance()->Ndim();
    long double rand_timings;

    Format::Format_A fmt(Nvol,Ndim);
    valarray<double> pj(fmt.size());

    FINE_TIMING_START(rand_timings);

    MPrand::mp_get_gauss(pj,rand);

    FINE_TIMING_END(rand_timings);
    _Message(TIMING_VERB_LEVEL, "[Timing] - MDutils::md_mom_su3 - RNG timing = "
	     << rand_timings << std::endl);

    pj *= sqrt(2.0); // to get Px,mu^a with the distribution exp(-1/2*P*P)
    pj /= 2.0;       // to get i*Px,mu =i*Px,mu^a*lambda^a/2 below

    
    for(int d=0; d< CommonPrms::instance()->Ndim(); ++d){
      for(int site=0; site<Nvol; ++site){
	
	double pjp[] = {0.0,   
			pj[fmt.index(2,site,d)]+sq3i*pj[fmt.index(7,site,d)],
			pj[fmt.index(1,site,d)], 
			pj[fmt.index(0,site,d)],          
			pj[fmt.index(4,site,d)], 
			pj[fmt.index(3,site,d)],          
			-pj[fmt.index(1,site,d)], 
			pj[fmt.index(0,site,d)],           
			0.0,  
			-pj[fmt.index(2,site,d)]+sq3i*pj[fmt.index(7,site,d)],
			pj[fmt.index(6,site,d)], 
			pj[fmt.index(5,site,d)],          
			-pj[fmt.index(4,site,d)], 
			pj[fmt.index(3,site,d)],           
			-pj[fmt.index(6,site,d)], 
			pj[fmt.index(5,site,d)],           
			0.0,  
			-2*sq3i*pj[fmt.index(7,site,d)]};
	
	P.data.set(P.format.cslice(0,site,d),valarray<double>(pjp,18));
      }
    }
  }

  void md_mom_suN(SUNRep<NC_>& Rep,GaugeField& P,const RandNum& rand){

    int Nvol = CommonPrms::instance()->Nvol();
    int Ndim = CommonPrms::instance()->Ndim();
    vector<int> gsite = SiteIndex::instance()->get_gsite();

    Format::Format_A fmt(Nvol,Ndim);
    valarray<double> pj(fmt.size());

    MPrand::mp_get_gauss(pj,rand,gsite,fmt);

    SUNmatrix<NC_> p, temp;    
    for(int d=0; d< Ndim; ++d){
      for(int site=0; site<Nvol; ++site){

	p.zero();
	for (int a = 0; a < ADJCOL; a++){
	  temp = Rep.lambda_fund(a).xI();
	  temp *= pj[fmt.index(a,site,d)];
	  p += temp;
	}

	P.data.set(P.format.cslice(0,site,d),p.getva());
      }
    }
  }

}
