/*!
  @file siteIndexEO.cpp
  @brief Declares the interface for site indexing EVEN-ODD
  Time-stamp: <2014-05-30 09:54:49 noaki>
*/
#include "siteIndex.hpp"
#include "Communicator/communicator.hpp"
#include "include/commonPrms.hpp"

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

int (SiteIndex::*SiteIndex::slice_dir[])(int, int)const = {&SiteIndex::slice_x,
							   &SiteIndex::slice_y,
							   &SiteIndex::slice_z,
							   &SiteIndex::slice_t};

int (SiteIndex::*SiteIndex::global_idx[])(int)const = {&SiteIndex::global_x,
						       &SiteIndex::global_y,
						       &SiteIndex::global_z,
						       &SiteIndex::global_t};


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

int SiteIndex::g_x(int gsite) const{ return gsite%Lx_;}
int SiteIndex::g_y(int gsite) const{ return (gsite/Lx_)%Ly_;}
int SiteIndex::g_z(int gsite) const{ return (gsite/LxLy_)%Lz_;}
int SiteIndex::g_t(int gsite) const{ return gsite/LxLyLz_;}

int SiteIndex::global_x(int x) const{
  return Communicator::instance()->ipe(XDIR)*Nx_+x;}

int SiteIndex::global_y(int y) const{ 
  return Communicator::instance()->ipe(YDIR)*Ny_+y;}

int SiteIndex::global_z(int z) const{ 
  return Communicator::instance()->ipe(ZDIR)*Nz_+z;}

int SiteIndex::global_t(int t) const{ 
  return Communicator::instance()->ipe(TDIR)*Nt_+t;}

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
  int n1 = (x&1)*((Ny_*Nz_*Nt_)/2+n)%((Ny_*Nz_*Nt_))+(1-x&1)*n;
  int parity = 2*n1/(Ny_*Nz_*Nt_);
  int idx = 2*n1-parity*(Ny_*Nz_*Nt_);
  return idx*Nxh_+ ((x+parity+idx%(Nz_*Ny_)/Ny_+idx%Ny_+idx/(Ny_*Nz_)) %2)*Nxh_+parity*Nvolh_+(x/2);}
int SiteIndex::slice_y(int y,int n) const{
  int n1 = (y&1)*((Nxh_*Nz_*Nt_)+n)%(Nx_*Nz_*Nt_)+(1-y&1)*n;
  return (n1/Nxh_)*Nxh_*Ny_+y*Nxh_+n1%Nxh_;}
int SiteIndex::slice_z(int z,int n) const{
  int n1 = (z&1)*((Nxh_*Ny_*Nt_)+n)%((Nx_*Ny_*Nt_))+(1-z&1)*n;
  return (n1/(Nxh_*Ny_))*Nxh_*Ny_*Nz_+z*Nxh_*Ny_+n1%(Nxh_*Ny_);}
int SiteIndex::slice_t(int t,int n) const{
  int n1 = (t&1)*((Nxh_*Ny_*Nz_)+n)%((Nx_*Ny_*Nz_))+(1-t&1)*n;
  return (n1/(Nxh_*Ny_*Nz_))*Nxh_*Ny_*Nz_*Nt_+n1%(Nxh_*Ny_*Nz_)+t*(Nxh_*Ny_*Nz_);}

int SiteIndex::slsize(int x,int dir)const{ return slsize_[dir];}

int SiteIndex::Bdir(int dir){return Bdir_[dir];}
const std::vector<int>& SiteIndex::get_gsite() const{ return global_site_;}
const std::vector<int>& SiteIndex::get_lsite() const{ return local_site_;}
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
  
  for(int site=0; site<Nvol_; ++site){
    local_site_.push_back(site);
    global_site_.push_back(c_x(site)+Nx_*nx 
			   +Lx*((c_y(site)+Ny_*ny)
				+Ly*((c_z(site)+Nz_*nz)
				     +Lz*(c_t(site)+Nt_*nt))));
  }
}
