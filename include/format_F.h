//----------------------------------------------------------------------
// format_F.h
//----------------------------------------------------------------------
#ifndef FORMAT_F_INCLUDED
#define FORMAT_F_INCLUDED

#ifndef SITEINDEX_INCLUDED
#include "Main/Geometry/siteIndex.h"
#endif

#include <valarray>

namespace Format{

  class Format_F{
  private:
    int Nvol_,Nex_,Nc_,Nd_;
    int Ndim_;
    int Nc2_;
    int Nin_;
  public:
    Format_F(int Nvol, int Nex=1):Nvol_(Nvol),Nex_(Nex),
				  Nc_(CommonPrms::instance()->Nc()),
				  Nd_(CommonPrms::instance()->Nd()),
				  Ndim_(CommonPrms::instance()->Ndim()),
				  Nc2_(2*Nc_),Nin_(Nc2_*Nd_){}

    int Nin() const {return Nin_;}
    int Nvol() const {return Nvol_;}
    int Nex() const {return Nex_;}
    int size() const {return Nin_*Nvol_*Nex_;}
    
    // get indices
    int index(int i, int site, int ex=0) const {
      return i +Nin_*site;
    }
    int index_r(int c, int s, int site, int ex=0) const { 
      return 2*(c +Nc_*s) +Nin_*site; 
    }
    int index_i(int c, int s, int site, int ex=0) const { 
      return 1+2*(c +Nc_*s) +Nin_*site; 
    }
    // get slices
    std::slice cplx_slice(int c, int s, int site, int ex=0) const {
      return std::slice(index_r(c,s,site), 2,1);
    }
    std::slice islice(int site, int ex=0) const {
      return std::slice(index(0,site), Nin_,1);
    }
    std::slice cslice(int s,int site,int ex=0) const {
      return std::slice(index_r(0,s,site), Nc2_,1);
    }
    std::gslice sslice(int c,int site,int ex=0) const {
      std::valarray<size_t> vsz_(2);
      std::valarray<size_t> vstr_(2);
      vsz_[0] = Nd_; vsz_[1] = 2;
      vstr_[0] = Nc2_; vstr_[1] = 1;

      return std::gslice(index_r(c,0,site), vsz_, vstr_);
    }
    
    std::gslice cs_slice(int c,int s) const {
      std::valarray<size_t> vsz_(2);
      std::valarray<size_t> vstr_(2);
      vsz_[0] = Nvol_; vsz_[1] = 2;
      vstr_[0] = Nin_; vstr_[1] = 1;

      return std::gslice(index_r(c,s,0), vsz_, vstr_);
    }

    const std::vector<int> get_sub(const std::vector<int>& sv)const{
      std::vector<int> sub;
      for(int e=0;e<Nex_;++e)
	for(int v=0;v<sv.size();++v)
	  for(int i=0;i<Nin_;++i)
	    sub.push_back(i+Nin_*(sv[v] +sv.size()*e));
      return sub;
    }
  };  
}  
#endif
