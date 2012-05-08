/*!
  @file siteIndexEO.cpp

  @brief Declares the interface for site indexing 

  EVEN-ODD
*/
#include "siteIndex.hpp"
#include "Communicator/communicator.h"

//std::vector<int> SiteIndex::global_site_;

SiteIndex::SiteIndex():Nx_(CommonPrms::instance()->Nx()),
		       Ny_(CommonPrms::instance()->Ny()),
		       Nz_(CommonPrms::instance()->Nz()),
		       Nt_(CommonPrms::instance()->Nt()),
		       Nvol_(CommonPrms::instance()->Nvol()),
		       NxNy_(Nx_*Ny_),NxNyNz_(Nx_*Ny_*Nz_),
		       Nvolh_(Nvol_/2),
		       Nxh_(Nx_/2),
		       Lx_(CommonPrms::instance()->Lx()),
		       Ly_(CommonPrms::instance()->Ly()),
		       Lz_(CommonPrms::instance()->Lz()),
		       Lt_(CommonPrms::instance()->Lt()),
		       LxLy_(Lx_*Ly_),LxLyLz_(Lx_*Ly_*Lz_),
		       Lvolh_(LxLyLz_*Lt_/2),
		       Lxh_(Lx_/2)
{
  setup_all();
}

SiteIndex* SiteIndex::instance(){
  static SiteIndex site_index;
  return &site_index;
}


int SiteIndex::site(int x,int y,int z,int t) const{
  return ((x+y+z+t) & 1)*Nvolh_ + (x >> 1) + Nxh_*y + Nxh_*Ny_*z + Nxh_*Ny_*Nz_*t;
}

int SiteIndex::c_x(int site) const{
  return ((site%Nxh_) << 1) + ((c_y(site)+c_z(site)+c_t(site) + (site/Nvolh_)) & 1);
}
int SiteIndex::c_y(int site) const{ return (site%(Nxh_*Ny_))/Nxh_;}
int SiteIndex::c_z(int site) const{ return (site%(Nxh_*Ny_*Nz_))/(Ny_*Nxh_);}
int SiteIndex::c_t(int site) const{ return (site%(Nxh_*Ny_*Nz_*Nt_))/(Nxh_*Ny_*Nz_);}


int SiteIndex::g_x(int gsite) const{
 return ((gsite%Lxh_) << 1) + ((g_y(gsite)+g_z(gsite)+g_t(gsite) + (gsite/Lvolh_)) & 1);
}
int SiteIndex::g_y(int gsite) const{ return (gsite%(Lxh_*Ly_))/Lxh_;}
int SiteIndex::g_z(int gsite) const{ return (gsite%(Lxh_*Ly_*Lz_))/(Ly_*Lxh_);}
int SiteIndex::g_t(int gsite) const{ return (gsite%(Lxh_*Ly_*Lz_*Lt_))/(Lxh_*Ly_*Lz_);}


int SiteIndex::global_x(int x){
  Communicator::instance()->ipe(XDIR)*Nx_+x;}

int SiteIndex::global_y(int y){ 
  Communicator::instance()->ipe(YDIR)*Ny_+y;}

int SiteIndex::global_z(int z){ 
  Communicator::instance()->ipe(ZDIR)*Nz_+z;}

int SiteIndex::global_t(int t){ 
  Communicator::instance()->ipe(TDIR)*Nt_+t;}

// indices with a step forward/backward (for bulk sites)
int SiteIndex::p_x(int site) const{ return (site+Nvolh_)%Nvol_ + (c_x(site) & 1);}
int SiteIndex::p_y(int site) const{ return (site+Nvolh_)%Nvol_+Nxh_;}
int SiteIndex::p_z(int site) const{ return (site+Nvolh_)%Nvol_+Nxh_*Ny_;}
int SiteIndex::p_t(int site) const{ return (site+Nvolh_)%Nvol_+Nxh_*Ny_*Nz_;}

int SiteIndex::m_x(int site) const{ return (site+Nvolh_)%Nvol_ - !(c_x(site) & 1);}
int SiteIndex::m_y(int site) const{ return (site+Nvolh_)%Nvol_-Nxh_;}
int SiteIndex::m_z(int site) const{ return (site+Nvolh_)%Nvol_-Nxh_*Ny_;}
int SiteIndex::m_t(int site) const{ return (site+Nvolh_)%Nvol_-Nxh_*Ny_*Nz_;}

int SiteIndex::slice_x(int x,int n) const{
  int parity = 2*n/(Ny_*Nz_*Nt_);
  int idx = 2*n-parity*(Ny_*Nz_*Nt_);
  return idx*Nxh_+ ((x+parity+idx%(Nz_*Ny_)/Nz_+idx%Ny_+idx/(Ny_*Nz_)) & 1)*Nxh_+parity*Nvolh_+(x >> 1);}
int SiteIndex::slice_y(int y,int n) const{
  return (n/Nxh_)*Nxh_*Ny_+y*Nxh_+n%Nxh_;}
int SiteIndex::slice_z(int z,int n) const{ return (n/(Nxh_*Ny_))*Nxh_*Ny_*Nz_+z*Nxh_*Ny_+n%(Nxh_*Ny_);}
int SiteIndex::slice_t(int t,int n) const{ return t*NxNyNz_+n;}

int SiteIndex::slsize(int x,int dir)const{ return slsize_[dir];}

int SiteIndex::Bdir(int dir){return Bdir_[dir];}
const std::vector<int>& SiteIndex::get_gsite() const{ return global_site_;}
int SiteIndex::get_gsite(int site) const{ return global_site_[site];}

void SiteIndex::setup_all() {
    int Bdir[Ndim_max_]= {Nx_-1,Ny_-1,Nz_-1,Nt_-1,};
    for(int i=0; i< Ndim_max_;++i) Bdir_.push_back(Bdir[i]);

    slsize_.push_back(Ny_*Nz_*Nt_);//x
    slsize_.push_back(Nx_*Nz_*Nt_);//y
    slsize_.push_back(Nx_*Ny_*Nt_);//z
    slsize_.push_back(Nx_*Ny_*Nz_);//t
    setup_global();
}

// setup of the global list vector
void SiteIndex::setup_global() {
  Communicator* commu = Communicator::instance();
  int nx = commu->ipe(XDIR);
  int ny = commu->ipe(YDIR);
  int nz = commu->ipe(ZDIR);
  int nt = commu->ipe(TDIR);

  int Lx = CommonPrms::instance()->Lx();
  int Ly = CommonPrms::instance()->Ly();
  int Lz = CommonPrms::instance()->Lz();
  
  for(int site=0; site<Nvol_; ++site)
    global_site_.push_back(c_x(site)+Nx_*nx 
			   +Lx*((c_y(site)+Ny_*ny)
				+Ly*((c_z(site)+Nz_*nz)
				     +Lz*(c_t(site)+Nt_*nt))));
}

