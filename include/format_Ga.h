//---------------------------------------------------------------------
// format_Ga.h
//---------------------------------------------------------------------
#ifndef FORMAT_GA_INCLUDED
#define FORMAT_GA_INCLUDED

#include "include/macros.hpp"
#include <vector>
#include <valarray>

namespace Format{

  class Format_Ga{
  private:
    const int Nvol_;
    const int Na_,Ndim_;
    const int Nin_;
    int Nex_;
  public:
    Format_Ga(int Nvol, int Nex=0):Nvol_(Nvol),
				   Na_(NC_*NC_-1),Ndim_(NDIM_),
				   Nex_(Ndim_),
				   Nin_(Na_*Na_){
      if(Nex) Nex_= Nex;
    }
    static int Nin(){return (NC_*NC_-1)*(NC_*NC_-1);}
    int Nvol() const {return Nvol_;}
    int Nex()  const {return Nex_;}
    int size() const {return Nin_*Nvol_*Nex_;}

    // get indices 
    int index(int aa,int site,int dir=0) const {
      return aa +Nin_*(site +Nvol_*dir);
    }
    int index_r(int a1,int a2,int site,int dir=0) const { 
      return Na_*a1 +a2 +Nin_*(site +Nvol_*dir); 
    }
    // get slices 
    std::slice islice(int site, int dir=0) const {
      return std::slice(index(0,site,dir),Nin_,1);
    }
    std::slice cslice(int in,int site,int dir=0) const { // "in" is dummy
      return std::slice(index(0,site,dir),Nin_,1);
    }
    std::slice ex_slice(int ex) const {
      return std::slice(index(0,0,ex),Nin_*Nvol_,1);
    }
    const std::valarray<size_t> get_sub(const std::vector<int>& sv)const{
      std::valarray<size_t> sub(Nin_*sv.size()*Nex_);
      int j=0;
      for(int e=0;e<Nex_;++e)
        for(int v=0;v<sv.size();++v)
          for(int i=0;i<Nin_;++i) sub[j++] = index(i,sv[v],e);
      return sub;
    }
    const std::valarray<size_t> get_sub(const std::vector<int>& sv,int e)const{
      std::valarray<size_t> sub(Nin_*sv.size());
      int j=0;
      for(int v=0;v<sv.size();++v)
	for(int i=0;i<Nin_;++i) sub[j++] = index(i,sv[v],e);
      return sub;
    }
  };
}

#endif
