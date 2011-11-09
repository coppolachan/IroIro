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
#include <cassert>

namespace MPrand{

  template<typename FMT>
  void mp_get(std::valarray<double>& rn,
	      const RandNum& rand,const FMT& fmt){

    int NP = CommonPrms::instance()->NP();
    int Nvol = CommonPrms::instance()->Nvol();
    int Lvol = CommonPrms::instance()->Lvol();
    int Nin = fmt.Nin();
    int Nex = fmt.Nex();

    std::valarray<double> Rn(NP*rn.size());
    rand.get(Rn);	

    assert(rn.size()==fmt.size());
    FMT Fmt(Lvol,Nex);
    
    for(int ex=0; ex<Nex;++ex){
      for(int site=0; site<Nvol; ++site){
	int gsite = SiteIndex::instance()->gsite(site);
	for(int in=0; in<Nin;++in){
	  rn[fmt.index(in,site,ex)] = Rn[Fmt.index(in,gsite,ex)];
	}
      }
    }    
  }
    
  template<typename FMT>
  void mp_get_gauss(std::valarray<double>& rn,
		    const RandNum& rand,const FMT& fmt){

    int NP = CommonPrms::instance()->NP();
    int Nvol = CommonPrms::instance()->Nvol();
    int Lvol = CommonPrms::instance()->Lvol();
    int Nin = fmt.Nin();
    int Nex = fmt.Nex();

    std::valarray<double> Rn(NP*rn.size());
    rand.get_gauss(Rn);	

    assert(rn.size()==fmt.size());
    FMT Fmt(Lvol,Nex);

    for(int ex=0; ex<fmt.Nex();++ex){
      for(int site=0; site<Nvol; ++site){
	int gsite = SiteIndex::instance()->gsite(site);

	for(int in=0; in<Nin;++in){
	  rn[fmt.index(in,site,ex)] = Rn[Fmt.index(in,gsite,ex)];
	}
      }
    }
  }

  void mp_get_gauss(std::valarray<double>& rn,const RandNum& rand);
  void mp_get(std::valarray<double>& rn,const RandNum& rand);
}

#endif 
