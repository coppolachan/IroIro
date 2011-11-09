//---------------------------------------------------------------------
// siteIndex_eo.cpp
//---------------------------------------------------------------------
#include "siteIndex_eo.h"
#include "Communicator/communicator.h"
#include <cassert>

using namespace std;

vector<vector<int> > SiteIndex_eo::ebdup_;
vector<vector<int> > SiteIndex_eo::ebdlw_;
vector<vector<int> > SiteIndex_eo::obdup_;
vector<vector<int> > SiteIndex_eo::obdlw_;
vector<vector<int> > SiteIndex_eo::ebprj_up_;
vector<vector<int> > SiteIndex_eo::ebprj_lw_;
vector<vector<int> > SiteIndex_eo::obprj_up_;
vector<vector<int> > SiteIndex_eo::obprj_lw_;

vector<int> SiteIndex_eo::esec_;
vector<int> SiteIndex_eo::osec_;
vector<int> SiteIndex_eo::ev_;
vector<int> SiteIndex_eo::ov_;

SiteIndex_eo* SiteIndex_eo::instance(){
  static SiteIndex_eo site_index;
  return &site_index;
}

void SiteIndex_eo::setup_eo(){

  int Nvol= CommonPrms::Nvol();
  int Ndim= CommonPrms::Ndim();
  
  int bsum = Communicator::instance()->ipe(0)*Nx_
            +Communicator::instance()->ipe(1)*Ny_
            +Communicator::instance()->ipe(2)*Nz_
            +Communicator::instance()->ipe(3)*Nt_;  
  int e=0;
  int o=0;

  ev_.resize(Nvol);
  ov_.resize(Nvol);

  for(int site=0;site<Nvol;++site){
    int eo =(bsum 
	     +idx_->x(site)
	     +idx_->y(site)
	     +idx_->z(site)
	     +idx_->t(site))%2;
    if(eo){ 
      osec_.push_back(site);
      ov_[site] = o;
      ++o;
    }else{
      esec_.push_back(site);
      ev_[site] = e;
      ++e;
    }
  }

  assert(esec_.size()==Nvol/2);
  assert(osec_.size()==Nvol/2);

  // boundary for the even sector
  ebdup_.resize(Ndim);
  ebdlw_.resize(Ndim);
  ebprj_up_.resize(Ndim);
  ebprj_lw_.resize(Ndim);

  for(int d=0;d<Ndim;++d){
    int up =0;
    int lw =0;
    ebprj_up_[d].resize(Nvol/2);
    ebprj_lw_[d].resize(Nvol/2);

    for(int hs=0;hs<Nvol/2;++hs){
      int cmp = idx_->cmp(esec_[hs],d);
      
      if(cmp==idx_->Bdir(d)){
	ebdup_[d].push_back(hs);
	ebprj_up_[d][hs] = up;
	++up;
      }
      if(cmp==0){
	ebdlw_[d].push_back(hs);
	ebprj_lw_[d][hs] = lw;
	++lw;
      }
    }
  }

  // boundary for the odd sector
  obdup_.resize(Ndim);
  obdlw_.resize(Ndim);
  obprj_up_.resize(Ndim);
  obprj_lw_.resize(Ndim);

  for(int d=0;d<Ndim;++d){
    int up =0;
    int lw =0;
    obprj_up_[d].resize(Nvol/2);
    obprj_lw_[d].resize(Nvol/2);

    for(int hs=0;hs<Nvol/2;++hs){
      int cmp = idx_->cmp(osec_[hs],d);
      
      if(cmp==idx_->Bdir(d)){
	obdup_[d].push_back(hs);
	obprj_up_[d][hs] = up;
	++up;
      }
      if(cmp==0){
	obdlw_[d].push_back(hs);
	obprj_lw_[d][hs] = lw;
	++lw;
      }
    }
  }
}    

int SiteIndex_eo::ex_p(int hs,int d){return ov_[idx_->x_p(esec_[hs],d)];}
int SiteIndex_eo::ex_m(int hs,int d){return ov_[idx_->x_m(esec_[hs],d)];}

int SiteIndex_eo::ox_p(int hs,int d){return ev_[idx_->x_p(osec_[hs],d)];}
int SiteIndex_eo::ox_m(int hs,int d){return ev_[idx_->x_m(osec_[hs],d)];}

const vector<int> SiteIndex_eo::ebdup(int d){return ebdup_[d];}  
const vector<int> SiteIndex_eo::ebdlw(int d){return ebdlw_[d];}  
const vector<int> SiteIndex_eo::obdup(int d){return obdup_[d];}  
const vector<int> SiteIndex_eo::obdlw(int d){return obdlw_[d];}  

int SiteIndex_eo::ebid_up(int hs,int d){return ebprj_up_[d][hs];}
int SiteIndex_eo::ebid_lw(int hs,int d){return ebprj_lw_[d][hs];}
int SiteIndex_eo::obid_up(int hs,int d){return obprj_up_[d][hs];}
int SiteIndex_eo::obid_lw(int hs,int d){return obprj_lw_[d][hs];}

int SiteIndex_eo::e_cmp(int hs,int d){return idx_->cmp(esec_[hs],d);}
int SiteIndex_eo::o_cmp(int hs,int d){return idx_->cmp(osec_[hs],d);}

