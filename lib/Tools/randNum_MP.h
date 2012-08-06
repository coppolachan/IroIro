//-----------------------------------------------------------------------------
// randNum_MP.h
//----------------------------------------------------------------------------
#ifndef RANDNUM_MP_INCLUDED
#define RANDNUM_MP_INCLUDED

#ifndef RANDNUM_INCLUDED
#include "randNum.h"
#endif

#include "include/commonPrms.h"
#include "Communicator/communicator.h"
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

    int Nvol = fmt.Nvol();
    int Nin = fmt.Nin();
    int Nex = fmt.Nex();

    int NP = CommonPrms::instance()->NP();
    int Lvol = NP*Nvol;

    std::valarray<double> Rn(NP*rn.size());
    rand.get(Rn);	

    FMT Fmt(Lvol,Nex);
    rn = Rn[Fmt.get_sub(gsite)];
  }

  template<typename FMT>
  void mp_get_gauss(std::valarray<double>& rn,
		    const RandNum& rand,
		    const std::vector<int>& gsite,
		    const FMT& fmt){
    rn=0.0;

    assert(rn.size()==fmt.size());

    int Nvol = fmt.Nvol();
    int Nin = fmt.Nin();
    int Nex = fmt.Nex();

    int NP = CommonPrms::instance()->NP();
    int Lvol = NP*Nvol;

    //CCIO::cout << "MP_GET_GAUSS valarray allocation\n";
      /*
    std::valarray<double> Rn(NP*rn.size());// too big.
    CCIO::cout << "MP_GET_GAUSS get_gauss\n";
    rand.get_gauss(Rn);	
    
    FMT Fmt(Lvol,Nex);
    CCIO::cout << "MP_GET_GAUSS distribute\n";
    rn = Rn[Fmt.get_sub(gsite)];

      */
    //    FMT Fmt(Lvol,Nex);
    std::valarray<double> Rn_source(rn.size());
    for (int node = 0; node < NP; ++node) {
      
      //CCIO::cout << "MP_GET_GAUSS get_gauss "<< node << "\n";
      if (Communicator::instance()->primaryNode()) {
	rand.get_gauss(Rn_source);	
      } 
      Communicator::instance()->sync();
      Communicator::instance()->send_1to1(rn, Rn_source, rn.size(), node, 0, node);
      //rn = Rn_node[Fmt.get_sub(gsite)];
      
    }
    
  }

  void mp_get(std::valarray<double>& rn,const RandNum& rand);
  void mp_get_gauss(std::valarray<double>& rn,const RandNum& rand);

}

#endif 
