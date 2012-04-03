//---------------------------------------------------------------------
// siteIndex.hpp
//---------------------------------------------------------------------
#ifndef SITEINDEX_INCLUDED
#define SITEINDEX_INCLUDED

#include "include/commonPrms.h"
#include "include/macros.hpp"
#include<iostream>
#include<vector>

enum site_dir{XDIR,YDIR,ZDIR,TDIR};

class SiteIndex{
private:
  enum{Ndim_max_ = NDIM_};

  int Nx_,Ny_,Nz_,Nt_;
  int NxNy_,NxNyNz_;
  int Nvol_;
  int Lx_,Ly_,Lz_,Lt_;
  int LxLy_,LxLyLz_;
  
  std::vector<int> Bdir_;
  std::vector<int> slsize_;
  static std::vector<int> global_site_;

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
	      LxLy_(Lx_*Ly_),LxLyLz_(Lx_*Ly_*Lz_){
    //
    int Bdir[Ndim_max_]= {Nx_-1,Ny_-1,Nz_-1,Nt_-1,};
    for(int i=0; i< Ndim_max_;++i) Bdir_.push_back(Bdir[i]);

    slsize_.push_back(Ny_*Nz_*Nt_);
    slsize_.push_back(Nx_*Nz_*Nt_);
    slsize_.push_back(Nx_*Ny_*Nt_);
    slsize_.push_back(Nx_*Ny_*Nz_);
    setup_global();
  }

  SiteIndex(const SiteIndex&){}
  SiteIndex& operator=(const SiteIndex&);
  void setup_lists();
  void setup_global();

public:
  static SiteIndex* instance();

  int site(int x,int y,int z,int t) const{
    return x +Nx_*y +NxNy_*z +NxNyNz_*t; }

  // components
  int c_x(int site) const{ return site%Nx_;}
  int c_y(int site) const{ return (site/Nx_)%Ny_;}
  int c_z(int site) const{ return (site/NxNy_)%Nz_;}
  int c_t(int site) const{ return site/NxNyNz_;}

  int g_x(int gsite) const{ return gsite%Lx_;}
  int g_y(int gsite) const{ return (gsite/Lx_)%Ly_;}
  int g_z(int gsite) const{ return (gsite/LxLy_)%Lz_;}
  int g_t(int gsite) const{ return gsite/LxLyLz_;}

  // indices with a step forward/backward (for bulk sites)
  int p_x(int site) const{ return site+1;}
  int p_y(int site) const{ return site+Nx_;}
  int p_z(int site) const{ return site+NxNy_;}
  int p_t(int site) const{ return site+NxNyNz_;}

  int m_x(int site) const{ return site-1;}
  int m_y(int site) const{ return site-Nx_;}
  int m_z(int site) const{ return site-NxNy_;}
  int m_t(int site) const{ return site-NxNyNz_;}

  int slice_x(int x,int n) const{ return n*Nx_+x;}
  int slice_y(int y,int n) const{ return (n/Nx_)*NxNy_+y*Nx_+n%Nx_;}
  int slice_z(int z,int n) const{ return (n/NxNy_)*NxNyNz_+z*NxNy_+n%NxNy_;}
  int slice_t(int t,int n) const{ return t*NxNyNz_+n;}

  int slsize(int x,int dir)const{ return slsize_[dir];}

  int Bdir(int dir){return Bdir_[dir];}
  const std::vector<int>& get_gsite() const{ return global_site_;}
  int get_gsite(int site) const{ return global_site_[site];}
};

#endif
