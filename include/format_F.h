//----------------------------------------------------------------------
/*! @file format_F.h
  @brief Defines the Format_Fermion classes
*/
//----------------------------------------------------------------------
#ifndef FORMAT_F_INCLUDED
#define FORMAT_F_INCLUDED

#include "Communicator/comm_io.hpp"
#include "include/macros.hpp"
#include <valarray>

namespace Format{

  template<int NCOL=NC_>
  class Format_Fermion{
  private:
    int Nvol_,Nex_,Nd_;
    int Ndim_;
    int Nin_;
  public:
    Format_Fermion(int Nvol, int Nex=1):Nvol_(Nvol),Nex_(Nex),
					Nd_(ND_),
					Ndim_(NDIM_),
					Nin_(2*NCOL*ND_){}
    static int Nin(){return 2*NCOL*ND_;}
    int Nvol() const {return Nvol_;}
    int Nex() const {return Nex_;}
    int size() const {return Nin_*Nvol_*Nex_;}
    
    // get indices
    inline int index(int i, int site, int ex=0) const {
      return i +Nin_*(site +Nvol_*ex);
    }

    inline int index_r(int c, int s, int site, int ex=0) const { 
      return 2*(c +NCOL*s) +Nin_*(site +Nvol_*ex); 
    }

    inline int index_i(int c, int s, int site, int ex=0) const { 
      return 1+2*(c +NCOL*s) +Nin_*(site +Nvol_*ex); 
    }
    // get slices
    std::slice cplx_slice(int c, int s, int site, int ex=0) const {
      return std::slice(index_r(c,s,site,ex), 2,1);
    }
    std::slice islice(int site, int ex=0) const {
      return std::slice(index(0,site,ex), Nin_,1);
    }
    std::slice cslice(int s,int site,int ex=0) const {
      return std::slice(index_r(0,s,site,ex), 2*NCOL,1);
    }
    std::slice ex_slice(int ex) const{
      return std::slice(index(0,0,ex),Nin_*Nvol_,1);
    }
    std::gslice sslice(int c,int site,int ex=0) const {
      std::valarray<size_t> vsz(2);
      std::valarray<size_t> vstr(2);
      vsz[0] = Nd_; vsz[1] = 2;
      vstr[0] = 2*NCOL; vstr[1] = 1;

      return std::gslice(index_r(c,0,site,ex),vsz,vstr);
    }
    
    std::gslice cs_slice(int c,int s,int ex=0) const {
      std::valarray<size_t> vsz(2);
      std::valarray<size_t> vstr(2);
      vsz[0] = Nvol_; vsz[1] = 2;
      vstr[0] = Nin_; vstr[1] = 1;

      return std::gslice(index_r(c,s,0,ex),vsz,vstr);
    }
    /*
    std::slice csr_slice(int c,int s,int ex=0) const{
      return std::slice(index_r(c,s,0,ex),Nin_,Nvol_);
    }
    std::slice csi_slice(int c,int s,int ex=0) const{
      return std::slice(index_i(c,s,0,ex),Nin_,Nvol_);
    }
    */
    /*! @brief returns subset of indices corresponding to given vector*/
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

  typedef Format_Fermion<NC_> Format_F;
  typedef Format_Fermion<NADJ_> Format_Fa;
}
  
#endif
