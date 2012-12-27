//---------------------------------------------------------------------
// format_A.h
//---------------------------------------------------------------------
#ifndef FORMAT_A_INCLUDED
#define FORMAT_A_INCLUDED

#include "macros.hpp"
#include <valarray>

namespace Format{

  class Format_A{
  private:
    int Nin_,Nvol_,Nex_;
  public:
    Format_A(int Nvol,int Nex=1)
      :Nvol_(Nvol),Nex_(Nex),Nin_(NC_*NC_-1){}

    static int Nin(){return NC_*NC_-1;}
    int Nvol() const {return Nvol_;}
    int Nex() const {return Nex_;}
    int size() const {return Nin_*Nvol_*Nex_;}

    // get indices 
    int index(int a, int site, int ex=0) const {
      return a +Nin_*(site +Nvol_*ex);
    }
    // get slices 
    std::slice islice(int site, int ex=0) const {
      return std::slice(index(0,site,ex),Nin_,1);
    }
    std::slice dir_slice(int ex) const {
      return std::slice(index(0,0,ex),Nin_*Nvol_,1);
    }

    const std::valarray<size_t> get_sub(const std::vector<int>& sv,int e)const{
      std::valarray<size_t> sub(Nin_*sv.size());
      int j=0; 
      for(int v=0;v<sv.size();++v)
	for(int i=0;i<Nin_;++i) sub[j++] = index(i,sv[v],e);
      return sub;
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
