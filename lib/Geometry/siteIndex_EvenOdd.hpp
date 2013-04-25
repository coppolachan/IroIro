//---------------------------------------------------------------------
// siteIndex_EvenOdd.hpp
//---------------------------------------------------------------------
#ifndef SITEINDEX_EVENODD_INCLUDED
#define SITEINDEX_EVENODD_INCLUDED

#include "include/commonPrms.hpp"
#include <cassert>
#include <vector>

class SiteIndex_EvenOdd{
private:
  enum{Ndim_max_ = NDIM_};

  int Nx_,Ny_,Nz_,Nt_;
  int Nvolh_;
  int NxNyNzh_,NxNyNth_,NxNzNth_,NyNzNth_;
  int NxNyh_,NxNzh_,NxNth_,NyNzh_;
  int Nxh_;
  int NyNz_;
  
  std::vector<int> Bdir_;
  static std::vector<int> esec_;
  static std::vector<int> osec_;
  static std::vector<int> global_even_;

  typedef int (SiteIndex_EvenOdd::*SliceSize)(int) const;
  SliceSize sls_xe;
  SliceSize sls_xo;

  int slsize_x_Veven(int)const{return NyNzNth_;}
  int slsize_xe_Vodd(int x)const{return (Ny_*Nz_*Nt_+1)/2-x%2;}
  int slsize_xo_Vodd(int x)const{return (Ny_*Nz_*Nt_-1)/2+x%2;}

  // constructor
  SiteIndex_EvenOdd():Nx_(CommonPrms::instance()->Nx()),
		      Ny_(CommonPrms::instance()->Ny()),
		      Nz_(CommonPrms::instance()->Nz()),
		      Nt_(CommonPrms::instance()->Nt()),
		      Nvolh_(CommonPrms::instance()->Nvol()/2),
		      NxNyNzh_(Nx_*Ny_*Nz_/2),NxNyNth_(Nx_*Ny_*Nt_/2),
		      NxNzNth_(Nx_*Nz_*Nt_/2),NyNzNth_(Ny_*Nz_*Nt_/2),
		      NxNyh_(Nx_*Ny_/2),NxNzh_(Nx_*Nz_/2),NxNth_(Nx_*Nt_/2),
		      NyNz_(Ny_*Nz_),NyNzh_(Ny_*Nz_/2),Nxh_(Nx_/2){
    assert(Nx_%2 ==0);
    int Bdir[Ndim_max_]= {Nx_-1,Ny_-1,Nz_-1,Nt_-1,};
    for(int i=0; i< Ndim_max_;++i) Bdir_.push_back(Bdir[i]);
    setup();
  }
  
  SiteIndex_EvenOdd(const SiteIndex_EvenOdd&){}
  SiteIndex_EvenOdd& operator=(const SiteIndex_EvenOdd&);
  void setup();

public:
  static SiteIndex_EvenOdd* instance(); 

  const std::vector<int>& esec()const{ return esec_;}
  const std::vector<int>& osec()const{ return osec_;}
  int esec(int hs)const{ return esec_[hs];}
  int osec(int hs)const{ return osec_[hs];}

  int c_xe(int hs)const{ return 2*(hs%Nxh_)+(c_t(hs)+c_z(hs)+c_y(hs))%2;}
  int c_xo(int hs)const{ return 2*(hs%Nxh_)+(c_t(hs)+c_z(hs)+c_y(hs)+1)%2;}
  int c_y(int hs)const{ return (hs%NxNyh_)/Nxh_;}
  int c_z(int hs)const{ return (hs%NxNyNzh_)/NxNyh_;}
  int c_t(int hs)const{ return hs/NxNyNzh_;}

  int p_xeo(int hs)const{ return hs +(c_t(hs)+c_z(hs)+c_y(hs)+1)%2;}
  int p_xoe(int hs)const{ return hs +(c_t(hs)+c_z(hs)+c_y(hs))%2;} 
  int p_y(int hs)const{ return hs +Nxh_; }
  int p_z(int hs)const{ return hs +NxNyh_; }
  int p_t(int hs)const{ return hs +NxNyNzh_; }

  int m_xeo(int hs)const{ return hs -(c_t(hs)+c_z(hs)+c_y(hs))%2; } 
  int m_xoe(int hs)const{ return hs -(c_t(hs)+c_z(hs)+c_y(hs)+1)%2;}
  int m_y(int hs)const{ return hs -Nxh_; }
  int m_z(int hs)const{ return hs -NxNyh_; }
  int m_t(int hs)const{ return hs -NxNyNzh_; }

  int slice_xe(int x,int n)const{ 
    return (2*n +(2*n/NyNz_+2*n%NyNz_/Ny_+2*n%Ny_+x)%2)*Nxh_+x/2;}	
  int slice_xo(int x,int n)const{ 
    return (2*n +(2*n/NyNz_+2*n%NyNz_/Ny_+2*n%Ny_+x+1)%2)*Nxh_+x/2;}	
  int slice_y(int y,int n)const{return (n/Nxh_)*NxNyh_+y*Nxh_+n%Nxh_;}
  int slice_z(int z,int n)const{return (n/NxNyh_)*NxNyNzh_+z*NxNyh_+n%NxNyh_;}
  int slice_t(int t,int n)const{return t*NxNyNzh_+n;}

  int slsize_xe(int x)const{return (this->*sls_xe)(x);}
  int slsize_xo(int x)const{return (this->*sls_xo)(x);}
  int slsize_y(int)const{return NxNzNth_;}
  int slsize_z(int)const{return NxNyNth_;}
  int slsize_t(int)const{return NxNyNzh_;}

  int Bdir(int dir){return Bdir_[dir];}
  const std::vector<int>& get_gsite()const{ return global_even_;}
  int get_gsite(int hs) const{return global_even_[hs];}

};

#endif
