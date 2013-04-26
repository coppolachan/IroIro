/*!--------------------------------------------------------------------------
 * @file antiPeriodicBC.hpp
 * @brief Definition of the AntiPeriodicBC class
 *-------------------------------------------------------------------------*/
#ifndef ANTIPERIODICBC_INCLUDED
#define ANTIPERIODICBC_INCLUDED

#include "include/pugi_interface.h"
#include "Tools/fieldUtils.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/sunAdjMatUtils.hpp"
#include "Geometry/siteMap.hpp"
#include <string.h>

template <typename F>
class AntiPeriodicBC{
  int dir_;
public:
  AntiPeriodicBC(int dir):dir_(dir){}
  void apply_bc(F& u)const;
  void apply_bc(F& ue,F& uo)const;
};

template <typename F>
void AntiPeriodicBC<F>::apply_bc(F& u)const{
  if(Communicator::ipe(dir_)==CommonPrms::NPE(dir_)-1){

    int Nbd = SiteIndex::instance()->Bdir(dir_);
    int slsize = SiteIndex::instance()->slsize(Nbd,dir_);

    for(int n=0; n<slsize; ++n){
      int site = SiteMap::shiftSite.xslice(Nbd,n,dir_);
      u.data.set(u.format.islice(site,dir_),
		 -FieldUtils::mat(u,site,dir_).getva());
    }
  }
}

template <typename F>
void AntiPeriodicBC<F>::apply_bc(F& ue,F& uo)const{
  if(Communicator::ipe(dir_)==CommonPrms::NPE(dir_)-1){

    int Nbd = SiteIndex_EvenOdd::instance()->Bdir(dir_);

    int slsize = SiteMap::shiftSite_eo.slice_size(Nbd,dir_);
    for(int n=0; n<slsize; ++n){
      int hs = SiteMap::shiftSite_eo.xslice(Nbd,n,dir_);
      ue.data.set(ue.format.islice(hs,dir_),
		  -FieldUtils::mat(ue,hs,dir_).getva());
    }
    slsize = SiteMap::shiftSite_oe.slice_size(Nbd,dir_);
    for(int n=0; n<slsize; ++n){
      int hs = SiteMap::shiftSite_oe.xslice(Nbd,n,dir_);
      uo.data.set(uo.format.islice(hs,dir_),
		  -FieldUtils::mat(uo,hs,dir_).getva());
    }
  }
}

/////////////
template<typename F>
AntiPeriodicBC<F>* createAPBC(XML::node bcnode){
  if(bcnode !=NULL){
    const char* bcname = bcnode.attribute("name").value();
    
    if(!strcmp(bcname,"periodic")) return NULL;
    if(!strcmp(bcname,"anti_periodic")){
      int dir;
      const char* dir_name = bcnode.attribute("dir").value();

      if(     !strcmp(dir_name,"X")) dir = XDIR;
      else if(!strcmp(dir_name,"Y")) dir = YDIR;
      else if(!strcmp(dir_name,"Z")) dir = ZDIR;
      else if(!strcmp(dir_name,"T")) dir = TDIR;
      else {
	CCIO::cout<<"No valid direction available\n";
	abort();
      }
      return new AntiPeriodicBC<F>(dir);
    }
    std::cerr<< "No valid name provided for boundary condition. Request by <"
	     << bcnode.name() << ">\n";
    abort();
  }
}
#endif
