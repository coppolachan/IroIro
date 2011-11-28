//---------------------------------------------------------------------
// siteIndex.cpp
//---------------------------------------------------------------------
#include <iostream>
#include "siteIndex.h"
#include "Communicator/communicator.h"

using namespace std;

vector<vector<int> > SiteIndex::bdup_;
vector<vector<int> > SiteIndex::bdlw_;

SiteIndex* SiteIndex::instance(){
  static SiteIndex site_index;
  return &site_index;
}

// arrays of the function pointers
int (SiteIndex::*SiteIndex::cmps[])(int site)const ={
  &SiteIndex::x, &SiteIndex::y, &SiteIndex::z, &SiteIndex::t,};

int (SiteIndex::*SiteIndex::xps[])(int site)const ={
  &SiteIndex::xp, &SiteIndex::yp, &SiteIndex::zp, &SiteIndex::tp,};

int (SiteIndex::*SiteIndex::xms[])(int site)const ={
  &SiteIndex::xm, &SiteIndex::ym, &SiteIndex::zm, &SiteIndex::tm,};

int (SiteIndex::*SiteIndex::xbds[])(int site)const ={
  &SiteIndex::xbd, &SiteIndex::ybd, &SiteIndex::zbd, &SiteIndex::tbd,};

void SiteIndex::setup_bdry(){
  int Nvol= CommonPrms::Nvol();
  int Ndim= CommonPrms::Ndim();

  bdup_.resize(Ndim);
  bdlw_.resize(Ndim);

  for(int d=0;d<Ndim;++d){
    bdup_[d].resize(Nvol/Ndir_[d]);
    bdlw_[d].resize(Nvol/Ndir_[d]);

    for(int site=0;site<Nvol;++site){
      if(cmp(site,d)==0)        bdlw_[d][x_b(site,d)]= site;
      if(cmp(site,d)==Bdir_[d]) bdup_[d][x_b(site,d)]= site;
    }
  }
}  

int SiteIndex::gsite(int site) const{
  
  CommonPrms* cprms = CommonPrms::instance();
  int Lx = cprms->Lx();
  int Ly = cprms->Ly();
  int Lz = cprms->Lz();
  int Lt = cprms->Lt();
  int Nx = cprms->Nx();
  int Ny = cprms->Ny();
  int Nz = cprms->Nz();
  int Nt = cprms->Nt();

  Communicator* commu = Communicator::instance();
  int nx = commu->ipe(0);
  int ny = commu->ipe(1);
  int nz = commu->ipe(2);
  int nt = commu->ipe(3);

  return x(site)+Nx*nx +Lx*((y(site)+Ny*ny)
			    +Ly*((z(site)+Nz*nz)
				 +Lz*(t(site)+Nt*nt)));
}

