//---------------------------------------------------------------------
// siteIndex.cpp
//---------------------------------------------------------------------
#include "siteIndex.h"
#include "Communicator/communicator.h"
#include <iostream>
#include <cassert>

using namespace std;

list_vec SiteIndex::bdup_;
list_vec SiteIndex::bdlw_;

list_vec SiteIndex::bdry_up_;
list_vec SiteIndex::bdry_lw_;
list_vec SiteIndex::bulk_up_;
list_vec SiteIndex::bulk_lw_;

vector<int> SiteIndex::gsite_;

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

  int Ndim= CommonPrms::Ndim();
  bdup_.resize(Ndim);
  bdlw_.resize(Ndim);
  bdry_up_.resize(Ndim);
  bdry_lw_.resize(Ndim);
  bulk_up_.resize(Ndim);
  bulk_lw_.resize(Ndim);

  for(int d=0;d<Ndim;++d){

    bdup_[d].resize(Nvol_/Ndir_[d]);
    for(int site=0;site<Nvol_;++site){
      if(cmp(site,d)==Bdir_[d]){
	bdup_[d][x_b(site,d)]= site;
	bdry_up_[d].push_back(site);
      }
      else bulk_up_[d].push_back(site);
    }
    bdlw_[d].resize(Nvol_/Ndir_[d]);
    for(int site=0;site<Nvol_;++site){
      if(cmp(site,d)==0){
        bdlw_[d][x_b(site,d)]= site;
	bdry_lw_[d].push_back(site);
      }
      else bulk_lw_[d].push_back(site);
    }
    assert(bdry_up_[d].size()==bdry_lw_[d].size());
  }
}  

void SiteIndex::setup_global() {
  Communicator* commu = Communicator::instance();
  int nx = commu->ipe(0);
  int ny = commu->ipe(1);
  int nz = commu->ipe(2);
  int nt = commu->ipe(3);

  for(int site=0; site<Nvol_; ++site)
    gsite_.push_back(x(site)+Nx_*nx 
		     +Lx_*((y(site)+Ny_*ny)
			   +Ly_*((z(site)+Nz_*nz)
				 +Lz_*(t(site)+Nt_*nt))));
}

