/*! 
  @file mdExec.cpp
  @brief utilities for MD including funcs to generate initial HMC momentum
 */
#include "mdExec.h"
#include "include/format_A.h"
#include "include/format_G.h"
#include "include/field.h"
#include "Tools/randNum_MP.h"

using namespace std;

namespace MDutils{
  void md_mom(Field& P,const RandNum& rand,const Format::Format_G& gf){
    if(CommonPrms::instance()->Nc() == 3){
      md_mom_su3(P,rand,gf);
    }else{
      throw "MD momentum P is not implemented for Nc_!=3.";
    }
  }

  void md_mom_su3(Field& P,const RandNum& rand,const Format::Format_G& gf){

    static const double sq3i = 1.0/sqrt(3.0);
    int Nvol = CommonPrms::instance()->Nvol();

    Format::Format_A fmt(Nvol);
    valarray<double> pj(fmt.size());

    vector<int> gsite;
    for(int site=0; site<Nvol; ++site)
      gsite.push_back(SiteIndex::instance()->gsite(site));

    MPrand::mp_get_gauss(pj,rand,gsite,fmt);
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
	
	P.set(gf.cslice(0,site,d),valarray<double>(pjp,18));
      }
    }
  }
}
