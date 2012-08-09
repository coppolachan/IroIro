//----------------------------------------------------------------------
/*! @file format_F.h
  @brief Defines the Format_F class
*/
//----------------------------------------------------------------------
#ifndef FORMAT_S_INCLUDED
#define FORMAT_S_INCLUDED

#include "Communicator/comm_io.hpp"
#include "include/macros.hpp"
#include <valarray>

namespace Format{

  class Format_S{
  private:
    int Nin_,Nvol_,Nex_;
    int Ndim_;
  public:
    Format_S(int Nvol,int Nex=1):Nvol_(Nvol),Nex_(Nex),
				 Ndim_(NDIM_),Nin_(2*NC_){}

    static int Nin(){return 2*NC_;}
    int Nvol() const {return Nvol_;}
    int Nex() const {return Nex_;}
    int size() const {return Nin_*Nvol_*Nex_;}
    
    // get indices
    inline int index(int i,int site,int ex=0) const{
      return i +Nin_*(site +Nvol_*ex);
    }
    inline int index_r(int c,int site,int ex=0) const{ 
      return 2*c +Nin_*(site +Nvol_*ex); 
    }
    inline int index_i(int c,int site,int ex=0) const{ 
      return 1+2*c +Nin_*(site +Nvol_*ex); 
    }
    // get slices
    std::slice cplx_slice(int c,int site,int ex=0) const{
      return std::slice(index_r(c,site,ex),2,1);
    }
    std::slice islice(int site,int ex=0) const{
      return std::slice(index(0,site,ex),Nin_,1);
    }
    std::slice cslice(int s,int site,int ex=0) const{
      return std::slice(index_r(0,site,ex),2*NC_,1);
    }
    std::slice ex_slice(int ex) const{
      return std::slice(index(0,0,ex),Nin_*Nvol_,1);
    }
    std::gslice c_slice(int c,int ex=0) const{
      return std::slice(index_r(c,0,ex),Nvol_,Nin_);
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
