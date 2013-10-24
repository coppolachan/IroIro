/*!
  @file siteIndex3d.hpp
  @brief Defines the interface for site indexing 
  At compile time is decided the if lexicographic or even-odd.
*/
#ifndef SITEINDEX3D_INCLUDED
#define SITEINDEX3D_INCLUDED

#include "include/commonPrms.hpp"
#include<vector>

class SiteIndex3d{
private:
  // Local dimensions
  int Nx_,Ny_,Nz_,NxNy_;
  
  std::vector<int> Bdir_;
  std::vector<int> slsize_;

  // hide copy constructors
  SiteIndex3d(const SiteIndex3d&){}  
  SiteIndex3d& operator=(const SiteIndex3d&);

public:
  SiteIndex3d():Nx_(CommonPrms::instance()->Nx()),
		Ny_(CommonPrms::instance()->Ny()),
		Nz_(CommonPrms::instance()->Nz()),
		NxNy_(Nx_*Ny_){
    slsize_.push_back(Ny_*Nz_); Bdir_.push_back(Nx_-1);
    slsize_.push_back(Nx_*Nz_); Bdir_.push_back(Ny_-1);
    slsize_.push_back(Nx_*Ny_); Bdir_.push_back(Nz_-1);
    slsize_.push_back(0); Bdir_.push_back(0);
  }
  int slice_x(int x,int n)const{return n*Nx_+x;}
  int slice_y(int y,int n)const{return (n/Nx_)*NxNy_+y*Nx_+n%Nx_;}
  int slice_z(int z,int n)const{return NxNy_*z +n;}
  int slice_t(int t,int n)const{return n;}

  int slsize(int dir)const{return slsize_[dir];}
  int Bdir(int dir) const{return Bdir_[dir];}
};

#endif
