//-----------------------------------------------------------------------------
// randNum_MP.h
//----------------------------------------------------------------------------
#ifndef RANDNUM_MP_INCLUDED
#define RANDNUM_MP_INCLUDED

#ifndef RANDNUM_INCLUDED
#include "randNum.h"
#endif

#ifndef SITEINDEX_INCLUDED
#include "Main/Geometry/siteIndex.h"
#endif

#include "Communicator/communicator.h"

#include <cassert>

namespace MPrand{

  template<typename FMT>
  void mp_get(std::valarray<double>& rn,
	      const RandNum& rand,
	      const std::vector<int>& gsite, 
	      const FMT& fmt){

    rn=0.0;
    assert(rn.size()==fmt.size());
    assert(rn.size()==gsite.size());

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
    assert(fmt.Nvol()==gsite.size());

    int Nvol = fmt.Nvol();
    int Nin = fmt.Nin();
    int Nex = fmt.Nex();

    int NP = CommonPrms::instance()->NP();
    int Lvol = NP*Nvol;

    std::valarray<double> Rn(NP*rn.size());
    rand.get_gauss(Rn);	
    
    FMT Fmt(Lvol,Nex);
    rn = Rn[Fmt.get_sub(gsite)];
  }

  void mp_get(std::valarray<double>& rn,const RandNum& rand);
  void mp_get_gauss(std::valarray<double>& rn,const RandNum& rand);

}

#endif 
