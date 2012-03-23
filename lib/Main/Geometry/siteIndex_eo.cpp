//---------------------------------------------------------------------
// siteIndex_eo.cpp
//---------------------------------------------------------------------
#include "siteIndex_eo.h"
#include "Communicator/communicator.h"
#include <cassert>

using namespace std;

vector<vector<int> > SiteIndex_eo::ebdup_;
vector<vector<int> > SiteIndex_eo::ebulk_up_;

vector<vector<int> > SiteIndex_eo::ebdlw_;
vector<vector<int> > SiteIndex_eo::ebulk_lw_;

vector<vector<int> > SiteIndex_eo::ebprj_up_;
vector<vector<int> > SiteIndex_eo::ebprj_lw_;

vector<int> SiteIndex_eo::esec_;
vector<int> SiteIndex_eo::ev_;

vector<vector<int> > SiteIndex_eo::obdup_;
vector<vector<int> > SiteIndex_eo::obulk_up_;

vector<vector<int> > SiteIndex_eo::obdlw_;
vector<vector<int> > SiteIndex_eo::obulk_lw_;

vector<vector<int> > SiteIndex_eo::obprj_up_;
vector<vector<int> > SiteIndex_eo::obprj_lw_;

vector<int> SiteIndex_eo::osec_;
vector<int> SiteIndex_eo::ov_;

vector<int> SiteIndex_eo::gev_;

SiteIndex_eo* SiteIndex_eo::instance(){
  static SiteIndex_eo site_index;
  return &site_index;
}

void SiteIndex_eo::setup_eo(){

  int bsum = Communicator::instance()->ipe(0)*Nx_
            +Communicator::instance()->ipe(1)*Ny_
            +Communicator::instance()->ipe(2)*Nz_
            +Communicator::instance()->ipe(3)*Nt_;  
  int e=0;
  int o=0;

  ev_.resize(Nvol_);
  ov_.resize(Nvol_);

  for(int site=0;site<Nvol_;++site){
    int eo =(bsum 
	     +idx_->x(site)
	     +idx_->y(site)
	     +idx_->z(site)
	     +idx_->t(site))%2;
    if(eo){ 
      osec_.push_back(site);
      ov_[site] = o++;
    }else{
      esec_.push_back(site);
      ev_[site] = e++;
    }
  }
  assert(esec_.size()==Nvol_/2);
  assert(osec_.size()==Nvol_/2);

  /*! global even-list vector */
  std::vector<int> Gev(Lvol_);
  e=0; 
  for(int gsite=0;gsite<Lvol_;++gsite){
    int eo =(idx_->gsx(gsite)
	    +idx_->gsy(gsite)
	    +idx_->gsz(gsite)
	    +idx_->gst(gsite))%2;
    
    if(!eo) Gev[gsite] = e++;
  }
  for(int hs=0; hs<Nvol_/2; ++hs) 
    gev_.push_back(Gev[idx_->get_gsite(esec_[hs])]);

  /*! boundary & bulk for the even sector */
  ebdup_.resize(Ndim_);      ebdlw_.resize(Ndim_);
  ebulk_up_.resize(Ndim_);   ebulk_lw_.resize(Ndim_);
  ebprj_up_.resize(Ndim_);   ebprj_lw_.resize(Ndim_);

  for(int d=0;d<Ndim_;++d){
    ebprj_up_[d].resize(Nvol_/2);
    int up =0;
    for(int hs=0;hs<Nvol_/2;++hs){
      if(idx_->cmp(esec_[hs],d) == idx_->Bdir(d)){
	ebdup_[d].push_back(hs);
	ebprj_up_[d][hs] = up++;
      }else{
	ebulk_up_[d].push_back(hs);
      }
    }
    ebprj_lw_[d].resize(Nvol_/2);
    int lw =0;
    for(int hs=0;hs<Nvol_/2;++hs){
      if(idx_->cmp(esec_[hs],d) == 0){
	ebdlw_[d].push_back(hs);
	ebprj_lw_[d][hs] = lw++;
      }else{
	ebulk_lw_[d].push_back(hs);
      }
    }
    assert(ebdup_[d].size()==ebdlw_[d].size());
  }
  /*! boundary & bulk for the odd sector */
  obdup_.resize(Ndim_);      obdlw_.resize(Ndim_);
  obulk_up_.resize(Ndim_);   obulk_lw_.resize(Ndim_);
  obprj_up_.resize(Ndim_);   obprj_lw_.resize(Ndim_);

  for(int d=0;d<Ndim_;++d){
    obprj_up_[d].resize(Nvol_/2);
    int up =0;
    for(int hs=0;hs<Nvol_/2;++hs){
      if(idx_->cmp(osec_[hs],d) == idx_->Bdir(d)){
	obdup_[d].push_back(hs);
	obprj_up_[d][hs] = up++;
      }else{
	obulk_up_[d].push_back(hs);
      }
    }
    obprj_lw_[d].resize(Nvol_/2);
    int lw =0;
    for(int hs=0;hs<Nvol_/2;++hs){
      if(idx_->cmp(osec_[hs],d) ==0){
	obdlw_[d].push_back(hs);
	obprj_lw_[d][hs] = lw++;
      }else{
	obulk_lw_[d].push_back(hs);
      }
    }
    assert(obdup_[d].size()==obdlw_[d].size());
  }
}    

int SiteIndex_eo::ex_p(int hs,int d)const{return ov_[idx_->x_p(esec_[hs],d)];}
int SiteIndex_eo::ex_m(int hs,int d)const{return ov_[idx_->x_m(esec_[hs],d)];}

int SiteIndex_eo::ox_p(int hs,int d)const{return ev_[idx_->x_p(osec_[hs],d)];}
int SiteIndex_eo::ox_m(int hs,int d)const{return ev_[idx_->x_m(osec_[hs],d)];}

const vector<int>& SiteIndex_eo::ebdup(int d)const{return ebdup_[d];}  
const vector<int>& SiteIndex_eo::ebdlw(int d)const{return ebdlw_[d];}  
const vector<int>& SiteIndex_eo::obdup(int d)const{return obdup_[d];}  
const vector<int>& SiteIndex_eo::obdlw(int d)const{return obdlw_[d];}  

int SiteIndex_eo::ebid_up(int hs,int d)const{return ebprj_up_[d][hs];}
int SiteIndex_eo::ebid_lw(int hs,int d)const{return ebprj_lw_[d][hs];}
int SiteIndex_eo::obid_up(int hs,int d)const{return obprj_up_[d][hs];}
int SiteIndex_eo::obid_lw(int hs,int d)const{return obprj_lw_[d][hs];}

int SiteIndex_eo::e_cmp(int hs,int d)const{return idx_->cmp(esec_[hs],d);}
int SiteIndex_eo::o_cmp(int hs,int d)const{return idx_->cmp(osec_[hs],d);}

