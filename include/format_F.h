//----------------------------------------------------------------------
/*! @file format_F.h

  @brief Defines the Format_F class

*/
//----------------------------------------------------------------------
#ifndef FORMAT_F_INCLUDED
#define FORMAT_F_INCLUDED

#include <valarray>
#include "Main/Geometry/siteIndex.h"
#include "Communicator/comm_io.hpp"

namespace Format{

  class Format_F{
  private:
    int Nvol_,Nex_,Nd_;
    int Ndim_;
    int Nin_;
  public:
    Format_F(int Nvol, int Nex=1):Nvol_(Nvol),Nex_(Nex),
				  Nd_(CommonPrms::instance()->Nd()),
				  Ndim_(CommonPrms::instance()->Ndim()),
				  Nin_(2*NC_*Nd_){}

    int Nin() const {return Nin_;}
    int Nvol() const {return Nvol_;}
    int Nex() const {return Nex_;}
    int size() const {return Nin_*Nvol_*Nex_;}
    
    // get indices
    inline int index(int i, int site, int ex=0) const {
      return i +Nin_*(site +Nvol_*ex);
    }

    inline int index_r(int c, int s, int site, int ex=0) const { 
      return 2*(c +NC_*s) +Nin_*(site +Nvol_*ex); 
    }

    inline int index_i(int c, int s, int site, int ex=0) const { 
      return 1+2*(c +NC_*s) +Nin_*(site +Nvol_*ex); 
    }
    // get slices
    std::slice cplx_slice(int c, int s, int site, int ex=0) const {
      return std::slice(index_r(c,s,site,ex), 2,1);
    }
    std::slice islice(int site, int ex=0) const {
      return std::slice(index(0,site,ex), Nin_,1);
    }
    std::slice cslice(int s,int site,int ex=0) const {
      return std::slice(index_r(0,s,site,ex), 2*NC_,1);
    }
    std::gslice sslice(int c,int site,int ex=0) const {
      std::valarray<size_t> vsz_(2);
      std::valarray<size_t> vstr_(2);
      vsz_[0] = Nd_; vsz_[1] = 2;
      vstr_[0] = 2*NC_; vstr_[1] = 1;

      return std::gslice(index_r(c,0,site,ex), vsz_, vstr_);
    }
    
    std::gslice cs_slice(int c,int s,int ex=0) const {
      std::valarray<size_t> vsz_(2);
      std::valarray<size_t> vstr_(2);
      vsz_[0] = Nvol_; vsz_[1] = 2;
      vstr_[0] = Nin_; vstr_[1] = 1;

      return std::gslice(index_r(c,s,0,ex), vsz_, vstr_);
    }

    /*! @brief returns subset of indices corresponding to given vector*/
    const std::valarray<size_t> get_sub(const std::vector<int>& sv)const{
      std::valarray<size_t> sub(Nin_*sv.size()*Nex_);
      int j=0;
      for(int e=0;e<Nex_;++e)
	for(int v=0;v<sv.size();++v)
	  for(int i=0;i<Nin_;++i) sub[j++] = index(i,sv[v],e);
      return sub;
    }
  };  
}
  
#endif
