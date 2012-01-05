//---------------------------------------------------------------------
// siteIndex.h
//---------------------------------------------------------------------
#ifndef SITEINDEX_INCLUDED
#define SITEINDEX_INCLUDED

#include "include/commonPrms.h"

#include<iostream>
#include<vector>

class SiteIndex{
private:
  enum{Ndim_max_ = 4};
  int Nx_,Ny_,Nz_,Nt_;
  int NxNy_,NxNyNz_;
  int Nvol_;
  int Lx_,Ly_,Lz_,Lt_;
  int LxLy_,LxLyLz_;
  std::vector<int> Bdir_;
  std::vector<int> Ndir_;
  std::vector<int> Vdir_;
  static std::vector<int> gsite_;
  static std::vector<std::vector<int> > bdup_;
  static std::vector<std::vector<int> > bdlw_;

  // indices with a step forward (for bulk sites)
  int xp(int site) const{ return site+1;}
  int yp(int site) const{ return site+Nx_;}
  int zp(int site) const{ return site+NxNy_;}
  int tp(int site) const{ return site+NxNyNz_;}

  // indices with a step backward (for bulk sites)
  int xm(int site) const{ return site-1;}
  int ym(int site) const{ return site-Nx_;}
  int zm(int site) const{ return site-NxNy_;}
  int tm(int site) const{ return site-NxNyNz_;}

  //indices on the boundaries
  int xbd(int site) const{ return site/Nx_;}
  int ybd(int site) const{ return site%Nx_+site/NxNy_*Nx_;}
  int zbd(int site) const{ return site%NxNy_+site/NxNyNz_*NxNy_;}
  int tbd(int site) const{ return site%NxNyNz_;}

  // constructor
  SiteIndex():Nx_(CommonPrms::instance()->Nx()),
	      Ny_(CommonPrms::instance()->Ny()),
	      Nz_(CommonPrms::instance()->Nz()),
	      Nt_(CommonPrms::instance()->Nt()),
	      Nvol_(CommonPrms::instance()->Nvol()),
	      NxNy_(Nx_*Ny_),NxNyNz_(Nx_*Ny_*Nz_),
	      Lx_(CommonPrms::instance()->Lx()),
	      Ly_(CommonPrms::instance()->Ly()),
	      Lz_(CommonPrms::instance()->Lz()),
	      Lt_(CommonPrms::instance()->Lt()),
	      LxLy_(Lx_*Ly_),LxLyLz_(Nx_*Ny_*Nz_){

    int Bdir[Ndim_max_]= {Nx_-1,Ny_-1,Nz_-1,Nt_-1,};
    int Ndir[Ndim_max_]= {Nx_,Ny_,Nz_,Nt_,};
    int Vdir[Ndim_max_]= {Nvol_/Nx_,Nvol_/Ny_,Nvol_/Nz_,Nvol_/Nt_,};

    for(int i=0; i< Ndim_max_;++i){
      Bdir_.push_back(Bdir[i]);
      Ndir_.push_back(Ndir[i]);
      Vdir_.push_back(Vdir[i]);
    }
    setup_bdry();
    setup_global();
  }

  SiteIndex(const SiteIndex&){}
  SiteIndex& operator=(const SiteIndex&){}
  void setup_bdry();
  void setup_global();
public:
  static SiteIndex* instance();
  // components
  int x(int site) const{ return site%Nx_;}
  int y(int site) const{ return (site/Nx_)%Ny_;}
  int z(int site) const{ return (site/NxNy_)%Nz_;}
  int t(int site) const{ return site/NxNyNz_;}

  int gx(int site) const{ return site%Lx_;}
  int gy(int site) const{ return (site/Lx_)%Ly_;}
  int gz(int site) const{ return (site/LxLy_)%Lz_;}
  int gt(int site) const{ return site/LxLyLz_;}

  int site(int x,int y,int z,int t) const{
    return x +Nx_*y +NxNy_*z +NxNyNz_*t; 
  }
  
  static int (SiteIndex::*cmps[])(int site) const;
  static int (SiteIndex::*xps[])(int site) const;
  static int (SiteIndex::*xms[])(int site) const;
  static int (SiteIndex::*xbds[])(int site) const;

  //// interfaces to handle the site index
  // component in each direction
  int cmp(int site,int dir)const{return (this->*SiteIndex::cmps[dir])(site);}

  // indices for one-step forward/backward in each direction
  int x_p(int site,int dir)const{return (this->*SiteIndex::xps[dir])(site);}
  int x_m(int site,int dir)const{return (this->*SiteIndex::xms[dir])(site);}

  // map from site to index on the boundary
  int x_b(int site,int dir)const{return (this->*SiteIndex::xbds[dir])(site);}

  int Bdir(int d)const{return Bdir_[d];}
  int Vdir(int d)const{return Vdir_[d];}
  const std::vector<int> bdup(int d)const{return bdup_[d];}
  const std::vector<int> bdlw(int d)const{return bdlw_[d];}

  const std::vector<int> get_gsite() const{ return gsite_;}
  int get_gsite(int site) const{ return gsite_[site];}
};

#endif
