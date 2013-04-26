//--------------------------------------------------------------
/*!@file dirac_staggeredLike.cpp
  @brief implementation of the utility functions for staggered fermions.
  Time-stamp: <2013-04-26 19:47:23 noaki>
*/
//-------------------------------------------------------------
#include "dirac_staggeredLike.hpp"
#include "Geometry/siteIndex.hpp"
#include "Geometry/siteIndex_EvenOdd.hpp"
using namespace std;

namespace Dstagg{

  void set_ksphase(valarray<double>& kse,valarray<double>& kso,int Nvol){
    // this function assumes kse and kso already have a correct size.
    kse = 1.0;
    kso = 1.0;
    for(int hs=0; hs<Nvol; ++hs){
      // setting for the even sector
      {
	int gs = SiteIndex::instance()->get_gsite(SiteIndex_EvenOdd::instance()->esec(hs));
	int x = SiteIndex::instance()->g_x(gs);
	int y = SiteIndex::instance()->g_y(gs);
	int z = SiteIndex::instance()->g_z(gs);
    
	kse[Nvol*YDIR +hs] *= 1.0-2.0*(x%2);
	kse[Nvol*ZDIR +hs] *= 1.0-2.0*((x+y)%2);
	kse[Nvol*TDIR +hs] *= 1.0-2.0*((x+y+z)%2);
      }
      // setting for the odd sector
      {
	int gs = SiteIndex::instance()->get_gsite(SiteIndex_EvenOdd::instance()->osec(hs)); 
	int x = SiteIndex::instance()->g_x(gs);
	int y = SiteIndex::instance()->g_y(gs);
	int z = SiteIndex::instance()->g_z(gs);
	
	kso[Nvol*YDIR +hs] *= 1.0-2.0*(x%2);
	kso[Nvol*ZDIR +hs] *= 1.0-2.0*((x+y)%2);
	kso[Nvol*TDIR +hs] *= 1.0-2.0*((x+y+z)%2);
      }
    }
  }
}
