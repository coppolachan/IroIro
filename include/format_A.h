//---------------------------------------------------------------------
// format_A.h
//---------------------------------------------------------------------
#ifndef FORMAT_A_INCLUDED
#define FORMAT_A_INCLUDED

#ifndef COMMONPRMS_INCLUDED
#include "commonPrms.h"
#endif

namespace Format{

  class Format_A{
  private:
    int Nvol_;
    int Nc_,Ndim_;
    int Nin_,Nex_;
  public:
    Format_A(int Nvol, int Nex=0):Nvol_(Nvol),
				  Nc_(CommonPrms::instance()->Nc()),
				  Ndim_(CommonPrms::instance()->Ndim()),
				  Nex_(Ndim_),
				  Nin_(Nc_*Nc_-1){
      if(Nex) Nex_= Nex;
    }
    int Nin() const {return Nin_;}
    int Nvol() const {return Nvol_;}
    int Nex() const {return Nex_;}
    int size() const {return Nin_*Nvol_*Nex_;}

    // get indices 
    int index(int a, int site, int dir=0) const {
      return a +Nin_*(site +Nvol_*dir);
    }
    // get slices 
    std::slice islice(int site, int dir=0) const {
      return std::slice(index(0,site,dir),Nin_,1);
    }
    std::slice dir_slice(int dir) const {
      return std::slice(index(0,0,dir),Nin_*Nvol_,1);
    }
  };
}

#endif
