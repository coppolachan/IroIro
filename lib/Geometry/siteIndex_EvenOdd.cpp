//---------------------------------------------------------------------
// siteIndex_EvenOdd.cpp
//---------------------------------------------------------------------
#include "siteIndex.hpp"
#include "siteIndex_EvenOdd.hpp"
#include "Communicator/communicator.hpp"
#include "Communicator/comm_io.hpp"

using namespace std;

vector<int> SiteIndex_EvenOdd::esec_;
vector<int> SiteIndex_EvenOdd::osec_;
vector<int> SiteIndex_EvenOdd::global_even_;

SiteIndex_EvenOdd* SiteIndex_EvenOdd::instance(){
  static SiteIndex_EvenOdd site_index;
  return &site_index;
}

// setup of the list vectors 
void SiteIndex_EvenOdd::setup(){

  if(Ny_*Nz_*Nt_%2){
    SiteIndex_EvenOdd::sls_xe = &SiteIndex_EvenOdd::slsize_xe_Vodd;
    SiteIndex_EvenOdd::sls_xo = &SiteIndex_EvenOdd::slsize_xo_Vodd;
  }else{
    SiteIndex_EvenOdd::sls_xe = &SiteIndex_EvenOdd::slsize_x_Veven;
    SiteIndex_EvenOdd::sls_xo = &SiteIndex_EvenOdd::slsize_x_Veven;
  }

  int bsum = Communicator::instance()->ipe(XDIR)*Nx_
            +Communicator::instance()->ipe(YDIR)*Ny_
            +Communicator::instance()->ipe(ZDIR)*Nz_
            +Communicator::instance()->ipe(TDIR)*Nt_;  
  
  SiteIndex *idx = SiteIndex::instance();

  for(int site=0;site<2*Nvolh_;++site){
    int eo =(bsum 
	     +idx->c_x(site) +idx->c_y(site)
	     +idx->c_z(site) +idx->c_t(site))%2;

    if(eo) osec_.push_back(site);
    else   esec_.push_back(site);
  }
  assert(esec_.size()==Nvolh_);
  assert(osec_.size()==Nvolh_);

  /*! global even-list vector */
  int Lvol = CommonPrms::instance()->Lvol();
  std::vector<int> Gev(Lvol);
  int e=0; 
  for(int gsite=0;gsite<Lvol;++gsite){
    int eo =( idx->g_x(gsite) +idx->g_y(gsite) 
	     +idx->g_z(gsite) +idx->g_t(gsite))%2;
    if(!eo) Gev[gsite] = e++;
    
  }
  for(int hs=0; hs<Nvolh_; ++hs) {
    global_even_.push_back(Gev[idx->get_gsite(esec_[hs])]);
  }
}    



