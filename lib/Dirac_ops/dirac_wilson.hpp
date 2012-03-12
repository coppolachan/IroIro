//----------------------------------------------------------------------
// dirac_wilson.hpp
//----------------------------------------------------------------------
#ifndef DIRAC_WILSON_INCLUDED
#define DIRAC_WILSON_INCLUDED

#include "include/format_F.h"
#include "include/format_G.h"
#include "include/pugi_interface.h"
#include "Main/Geometry/siteIndex_eo.h"

#include "Tools/sunMat.hpp"
#include "Tools/sunVec.hpp"
#include "dirac.hpp"

typedef Format::Format_F ffmt_t;
typedef Format::Format_G gfmt_t;

namespace Dw{
  struct EOtag{};    
  struct OEtag{};    

  double read_mass(const XML::node& node);
}

class Dirac_Wilson: public DiracWilsonLike{

private:
  const Field* const u_;
  double kpp_;
  int Nvol_;
  int Ndim_;
  int Nd_;
  //int boundary[Ndim];
  // this is temporary setting.

  const ffmt_t ff_;
  const gfmt_t gf_;

  const size_t fsize_;
  const size_t gsize_;

  const list_vec& bdry_mup_;
  const list_vec& bdry_mlw_;
  const list_vec& bdry_pup_;
  const list_vec& bdry_plw_;

  const list_vec& bulk_mup_;
  const list_vec& bulk_mlw_;
  const list_vec& bulk_pup_;
  const list_vec& bulk_plw_;

  void mult_xp(Field&,const Field&)const;
  void mult_yp(Field&,const Field&)const;
  void mult_zp(Field&,const Field&)const;
  void mult_tp(Field&,const Field&)const;

  void mult_xm(Field&,const Field&)const;
  void mult_ym(Field&,const Field&)const;
  void mult_zm(Field&,const Field&)const;
  void mult_tm(Field&,const Field&)const;

  static void(Dirac_Wilson::*mult_p[])(Field&,const Field&)const;
  static void(Dirac_Wilson::*mult_m[])(Field&,const Field&)const;

  inline int r0(int c)const{return 2*c;}
  inline int r1(int c)const{return 2*(NC_+c);}
  inline int r2(int c)const{return 2*(2*NC_+c);}
  inline int r3(int c)const{return 2*(3*NC_+c);} 

  inline int i0(int c)const{return 2*c+1;}
  inline int i1(int c)const{return 2*(NC_+c)+1;}
  inline int i2(int c)const{return 2*(2*NC_+c)+1;}
  inline int i3(int c)const{return 2*(3*NC_+c)+1;} 

  inline int re(int c1,int c2)const{return 2*(NC_*c1+c2);}
  inline int im(int c1,int c2)const{return 2*(NC_*c1+c2)+1;}

  void mult_full(Field&,const Field&)const;   /*! @brief  -kpp*D*f */
  void mult_offdiag(Field&,const Field&)const;/*! @brief  (1-kpp*D)*f */

  int gsite(int site)const {return site;}
  int esec(int site)const {return SiteIndex_eo::instance()->esec(site);}
  int osec(int site)const {return SiteIndex_eo::instance()->osec(site);}

  int xp_full(int site,int dir)const {
    return SiteIndex::instance()->x_p(site,dir);}
  int xp_eo(int hs,int dir)const {
    return SiteIndex_eo::instance()->ex_p(hs,dir);}
  int xp_oe(int hs,int dir)const {
    return SiteIndex_eo::instance()->ox_p(hs,dir);}

  int xm_full(int site,int dir)const {
    return SiteIndex::instance()->x_m(site,dir);}
  int xm_eo(int hs,int dir)const {
    return SiteIndex_eo::instance()->ex_m(hs,dir);}
  int xm_oe(int hs,int dir)const {
    return SiteIndex_eo::instance()->ox_m(hs,dir);}

  int(Dirac_Wilson::*gp)(int)const;
  int(Dirac_Wilson::*gm)(int)const;
  int(Dirac_Wilson::*x_p)(int,int)const;
  int(Dirac_Wilson::*x_m)(int,int)const;
  void(Dirac_Wilson::*mult_core)(Field&,const Field&)const;

  Dirac_Wilson(const Dirac_Wilson&); /*!< @brief simple copy is prohibited.*/

public:
  /*! @brief constructor to create instance with e/o site indexing */
  Dirac_Wilson(double mass,const Field* u,Dw::EOtag)
    :kpp_(0.5/(4.0+mass)),u_(u),
     Nvol_(CommonPrms::instance()->Nvol()/2),
     Ndim_(CommonPrms::instance()->Ndim()),
     Nd_(CommonPrms::instance()->Nd()),
     gp(&Dirac_Wilson::esec),gm(&Dirac_Wilson::osec),
     x_p(&Dirac_Wilson::xp_eo),x_m(&Dirac_Wilson::xm_eo),
     mult_core(&Dirac_Wilson::mult_offdiag),
     ff_(Nvol_),  fsize_(ff_.size()),
     gf_(2*Nvol_),gsize_(gf_.size()),
     bdry_mlw_(SiteIndex_eo::instance()->ebdry_lw()),
     bdry_plw_(SiteIndex_eo::instance()->obdry_lw()),
     bdry_mup_(SiteIndex_eo::instance()->obdry_up()),
     bdry_pup_(SiteIndex_eo::instance()->ebdry_up()),
     bulk_mlw_(SiteIndex_eo::instance()->ebulk_lw()),
     bulk_plw_(SiteIndex_eo::instance()->obulk_lw()),
     bulk_mup_(SiteIndex_eo::instance()->obulk_up()),
     bulk_pup_(SiteIndex_eo::instance()->ebulk_up()){}

  Dirac_Wilson(double mass,const Field* u,Dw::OEtag)
    :kpp_(0.5/(4.0+mass)),u_(u),
     Nvol_(CommonPrms::instance()->Nvol()/2),
     Ndim_(CommonPrms::instance()->Ndim()),
     Nd_(CommonPrms::instance()->Nd()),
     gp(&Dirac_Wilson::osec),gm(&Dirac_Wilson::esec),
     x_p(&Dirac_Wilson::xp_oe),x_m(&Dirac_Wilson::xm_oe),
     mult_core(&Dirac_Wilson::mult_offdiag),
     ff_(Nvol_), fsize_(ff_.size()),
     gf_(2*Nvol_),gsize_(gf_.size()),
     bdry_mlw_(SiteIndex_eo::instance()->obdry_lw()),
     bdry_plw_(SiteIndex_eo::instance()->ebdry_lw()),
     bdry_mup_(SiteIndex_eo::instance()->ebdry_up()),
     bdry_pup_(SiteIndex_eo::instance()->obdry_up()),
     bulk_mlw_(SiteIndex_eo::instance()->obulk_lw()),
     bulk_plw_(SiteIndex_eo::instance()->ebulk_lw()),
     bulk_mup_(SiteIndex_eo::instance()->ebulk_up()),
     bulk_pup_(SiteIndex_eo::instance()->obulk_up()){}

  /*! @brief constructor to create instance with normal site indexing
   */
  Dirac_Wilson(double mass,const Field* u)
    :kpp_(0.5/(4.0+mass)),u_(u),
     Nvol_(CommonPrms::instance()->Nvol()),
     Ndim_(CommonPrms::instance()->Ndim()),
     Nd_(CommonPrms::instance()->Nd()),
     gp(&Dirac_Wilson::gsite),gm(&Dirac_Wilson::gsite),
     x_p(&Dirac_Wilson::xp_full),x_m(&Dirac_Wilson::xm_full),
     mult_core(&Dirac_Wilson::mult_full),
     ff_(Nvol_),fsize_(ff_.size()),
     gf_(Nvol_),gsize_(gf_.size()),
     bdry_mlw_(SiteIndex::instance()->bdry_lw()),
     bdry_plw_(SiteIndex::instance()->bdry_lw()),
     bdry_mup_(SiteIndex::instance()->bdry_up()),
     bdry_pup_(SiteIndex::instance()->bdry_up()),
     bulk_mlw_(SiteIndex::instance()->bulk_lw()),
     bulk_plw_(SiteIndex::instance()->bulk_lw()),
     bulk_mup_(SiteIndex::instance()->bulk_up()),
     bulk_pup_(SiteIndex::instance()->bulk_up()){}

  Dirac_Wilson(const XML::node& node,const Field* u)
    :u_(u),
     Nvol_(CommonPrms::instance()->Nvol()),
     Ndim_(CommonPrms::instance()->Ndim()),
     Nd_(CommonPrms::instance()->Nd()),
     gp(&Dirac_Wilson::gsite),gm(&Dirac_Wilson::gsite),
     x_p(&Dirac_Wilson::xp_full),x_m(&Dirac_Wilson::xm_full),
     mult_core(&Dirac_Wilson::mult_full),
     ff_(Nvol_),fsize_(ff_.size()),
     gf_(Nvol_),gsize_(gf_.size()),
     bdry_mlw_(SiteIndex::instance()->bdry_lw()),
     bdry_plw_(SiteIndex::instance()->bdry_lw()),
     bdry_mup_(SiteIndex::instance()->bdry_up()),
     bdry_pup_(SiteIndex::instance()->bdry_up()),
     bulk_mlw_(SiteIndex::instance()->bulk_lw()),
     bulk_plw_(SiteIndex::instance()->bulk_lw()),
     bulk_mup_(SiteIndex::instance()->bulk_up()),
     bulk_pup_(SiteIndex::instance()->bulk_up()){
    //
    double mass;
    XML::read(node, "mass", mass);
    kpp_= 0.5/(4.0+mass);
  }

  virtual ~Dirac_Wilson(){}
  
  size_t fsize() const{return fsize_;}
  size_t gsize() const{return gsize_;}

  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;

  ////////////////////////////////////////Preconditioned versions
  // Wilson operator has no defined preconditioner now 
  const Field mult_prec     (const Field&f)const{return f;}
  const Field mult_dag_prec (const Field&f)const{return f;}
  const Field left_prec     (const Field&f)const{return f;}
  const Field right_prec    (const Field&f)const{return f;}
  const Field left_dag_prec (const Field&f)const{return f;}
  const Field right_dag_prec(const Field&f)const{return f;}
  //////////////////////////////////////////////////////////////

  const Field gamma5(const Field&) const;
  const Field proj_p(const Field&) const;
  const Field proj_m(const Field&) const;

  const Field md_force(const Field& , const Field&)const;
  void md_force_p(Field&,const Field&,const Field&)const;
  void md_force_m(Field&,const Field&,const Field&)const;
  
  double getKappa() const {return kpp_;}  
  const ffmt_t get_fermionFormat() const {return ff_;}
  const std::vector<int> get_gsite() const;

  void update_internal_state(){}
};
#endif
