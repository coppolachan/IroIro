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

  // even sector
  static list_vec ebdup_;
  static list_vec ebdlw_;
  static list_vec ebulk_up_;
  static list_vec ebulk_lw_;
  static list_vec ebprj_up_;
  static list_vec ebprj_lw_;
  static std::vector<int> esec_;
  static std::vector<int> ev_;

  // odd sector
  static list_vec obdup_;
  static list_vec obdlw_;
  static list_vec obulk_up_;
  static list_vec obulk_lw_;
  static list_vec obprj_up_;
  static list_vec obprj_lw_;
  static std::vector<int> osec_;
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
  int ex_p(int hs,int d)const;
  int ex_m(int hs,int d)const;
  // odd -> even
  int ox_p(int hs,int d)const;
  int ox_m(int hs,int d)const;

  // projection to the boundary space within the half field.
  const std::vector<int>& ebdup(int d) const;
  const std::vector<int>& ebdlw(int d) const;

  const std::vector<int>& obdup(int d) const;
  const std::vector<int>& obdlw(int d) const;
  
  // index in the boundary subspace
  int ebid_up(int hs,int d) const;
  int ebid_lw(int hs,int d) const;
  int obid_up(int hs,int d) const;
  int obid_lw(int hs,int d) const;

  int Ve_dir(int d)const{return ebdup_[d].size();}  
  int Vo_dir(int d)const{return obdup_[d].size();}  

  // component in the given direction d.
  int e_cmp(int hs,int d) const;
  int o_cmp(int hs,int d) const;

  int Bdir(int d)const {return idx_->Bdir(d);}

  const std::vector<int>& esec()const {return esec_;}
  const std::vector<int>& osec()const {return osec_;}
  int esec(int site)const {return esec_[site];}
  int osec(int site)const {return osec_[site];}

  const list_vec& ebdry_up()const {return ebdup_;}
  const list_vec& ebdry_lw()const {return ebdlw_;}
  const list_vec& ebulk_up()const {return ebulk_up_;}
  const list_vec& ebulk_lw()const {return ebulk_lw_;}

  const list_vec& obdry_up()const {return obdup_;}
  const list_vec& obdry_lw()const {return obdlw_;}
  const list_vec& obulk_up()const {return obulk_up_;}
  const list_vec& obulk_lw()const {return obulk_lw_;}

  const std::vector<int>& get_gsite() const{ return gev_;}
  int get_gsite(int hs) const{return gev_[hs];}
};

#endif
