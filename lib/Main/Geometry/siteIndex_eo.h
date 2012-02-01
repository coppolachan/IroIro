//---------------------------------------------------------------------
// siteIndex_eo.h
//---------------------------------------------------------------------
#ifndef SITEINDEX_EO_INCLUDED
#define SITEINDEX_EO_INCLUDED

#ifndef SITEINDEX_INCLUDED
#include "siteIndex.h"
#endif

class SiteIndex_eo{
  
private:
  SiteIndex* idx_;
  int Nx_,Ny_,Nz_,Nt_;
  int Nvol_,Ndim_,Lvol_;
  static std::vector<std::vector<int> > ebdup_;
  static std::vector<std::vector<int> > ebdlw_;
  static std::vector<std::vector<int> > obdup_;
  static std::vector<std::vector<int> > obdlw_;

  static std::vector<std::vector<int> > ebprj_up_;
  static std::vector<std::vector<int> > ebprj_lw_;
  static std::vector<std::vector<int> > obprj_up_;
  static std::vector<std::vector<int> > obprj_lw_;

  static std::vector<int> esec_;
  static std::vector<int> osec_;

  static std::vector<int> ev_;
  static std::vector<int> ov_;

  static std::vector<int> gev_;

  // constructor
  SiteIndex_eo():Nx_(CommonPrms::instance()->Nx()),
		 Ny_(CommonPrms::instance()->Ny()),
		 Nz_(CommonPrms::instance()->Nz()),
		 Nt_(CommonPrms::instance()->Nt()),
		 Nvol_(CommonPrms::instance()->Nvol()),
		 Lvol_(CommonPrms::instance()->Lvol()),
		 Ndim_(CommonPrms::instance()->Ndim()),
		 idx_(SiteIndex::instance()){ setup_eo();}
  
  SiteIndex_eo(const SiteIndex_eo&){}
  SiteIndex_eo& operator=(const SiteIndex_eo&);
  void setup_eo();

public:
  static SiteIndex_eo* instance();

  // index of one-step forward/backward in the given direction d.
  // even -> odd
  int ex_p(int hs,int d);
  int ex_m(int hs,int d);
  // odd -> even
  int ox_p(int hs,int d);
  int ox_m(int hs,int d);

  // projection to the boundary space within the half field.
  const std::vector<int> ebdup(int d);
  const std::vector<int> ebdlw(int d);

  const std::vector<int> obdup(int d);
  const std::vector<int> obdlw(int d);
  
  // index in the boundary subspace
  int ebid_up(int hs,int d);
  int ebid_lw(int hs,int d);
  int obid_up(int hs,int d);
  int obid_lw(int hs,int d);

  int Ve_dir(int d)const{return ebdup_[d].size();}  
  int Vo_dir(int d)const{return obdup_[d].size();}  

  // component in the given direction d.
  int e_cmp(int hs,int d);
  int o_cmp(int hs,int d);

  int Bdir(int d){return idx_->Bdir(d);}

  const std::vector<int>& esec(){return esec_;}
  const std::vector<int>& osec(){return osec_;}
  static int esec(int site){return esec_[site];}
  static int osec(int site){return osec_[site];}

  const std::vector<int> get_gsite() const{ return gev_;}
  int get_gsite(int hs) const{return gev_[hs];}
};

#endif
