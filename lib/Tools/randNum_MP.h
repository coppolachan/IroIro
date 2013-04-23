//-----------------------------------------------------------------------------
// randNum_MP.h
//----------------------------------------------------------------------------
#ifndef RANDNUM_MP_INCLUDED
#define RANDNUM_MP_INCLUDED

#include "randNum.h"

#include "include/commonPrms.h"
#include "Communicator/communicator.hpp"
#include "Communicator/comm_io.hpp"
#include <cassert>

namespace MPrand{

  template<typename FMT>
  void mp_get(std::valarray<double>& rn,
	      const RandNum& rand,
	      const std::vector<int>& gsite, 
	      const FMT& fmt){
    rn=0.0;
    assert(rn.size()==fmt.size());

    int NP = CommonPrms::instance()->NP();
    std::valarray<double> Rn(NP*rn.size());

    rand.get(Rn);	

    FMT Fmt(NP*fmt.Nvol(),fmt.Nex());
    rn = Rn[Fmt.get_sub(gsite)];
  }

  
  template<typename FMT>
  void mp_get_gauss(std::valarray<double>& rn,
		    const RandNum& rand,
		    const std::vector<int>& gsite,
		    const FMT& fmt){
    rn=0.0;

    assert(rn.size()==fmt.size());
    
    if (rand.parallel_safe()){
      rand.get_gauss(rn);
      Communicator::instance()->sync();	
    } else {
      
      int NP = CommonPrms::instance()->NP();
      
      std::valarray<double> Rn_source(rn.size());
      for(int node=0; node<NP; ++node) {
	if(Communicator::instance()->primaryNode()) 
	  rand.get_gauss(Rn_source);	
	
	Communicator::instance()->sync();
	Communicator::instance()->send_1to1(rn,Rn_source,rn.size(),node,0,node);
      }
    }
  }
  
  
  void mp_get(std::valarray<double>&,const RandNum&);
  void mp_get_gauss(std::valarray<double>&,const RandNum&);
  
  
}

#endif 
