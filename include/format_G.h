//---------------------------------------------------------------------
// format_G.h
//---------------------------------------------------------------------
#ifndef FORMAT_G_INCLUDED
#define FORMAT_G_INCLUDED

#include "include/macros.hpp"
#include <vector>
#include <valarray>

namespace Format{

  class Format_G{
  private:
    int Nvol_;
    int Nc_,Ndim_;
    int Nin_,Nex_;
  public:
    Format_G(int Nvol, int Nex=0):Nvol_(Nvol),
				  Nc_(NC_),Ndim_(NDIM_),
				  Nex_(Ndim_),
				  Nin_(2*NC_*NC_){
      if(Nex) Nex_= Nex;
    }
    static int Nin(){return 2*NC_*NC_;}
    int Nvol() const {return Nvol_;}
    int Nex()  const {return Nex_;}
    int size() const {return Nin_*Nvol_*Nex_;}

    // get indices 
    int index(int cc, int site, int dir=0) const {
      return cc +Nin_*(site +Nvol_*dir);
    }
    int index_r(int c1, int c2, int site, int dir=0) const { 
      return 2*(Nc_*c1 +c2) +Nin_*(site +Nvol_*dir); 
    }
    int index_i(int c1, int c2, int site, int dir=0) const { 
      return 1 +2*(Nc_*c1 +c2) +Nin_*(site +Nvol_*dir); 
    }
    // get slices 
    std::slice cplx_slice(int c1, int c2, int site, int dir=0) const {
      return std::slice(index_r(c1,c2,site,dir),2,1);
    }
    std::slice islice(int site, int dir=0) const {
      return std::slice(index(0,site,dir),Nin_,1);
    }
    std::slice cslice(int in,int site,int dir=0) const { // "in" is dummy
      return std::slice(index(0,site,dir),Nin_,1);
    }
    std::slice dir_slice(int dir) const {
      return std::slice(index(0,0,dir),Nin_*Nvol_,1);
    }
    std::gslice dslice(int site, int c1, int c2) const {
      std::valarray<std::size_t> vsz_(2);
      std::valarray<std::size_t> vstr_(2);
      vsz_[0] = Nex_;        vsz_[1] =2;
      vstr_[0] = Nin_*Nvol_; vstr_[1] =1;
      
      return std::gslice(index_r(c1,c2,site,0),vsz_,vstr_);
    }

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
