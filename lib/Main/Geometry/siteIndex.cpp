//---------------------------------------------------------------------
// siteIndex.cpp
//---------------------------------------------------------------------
#include "siteIndex.hpp"
#include "Communicator/communicator.h"
#include <iostream>
#include <cassert>

using namespace std;
vector<int> SiteIndex::global_site_;

SiteIndex* SiteIndex::instance(){
  static SiteIndex site_index;
  return &site_index;
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

