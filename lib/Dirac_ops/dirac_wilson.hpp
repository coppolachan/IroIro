/*! @file dirac_wilson.hpp
 * @brief Dirac_Wilson class 
 Time-stamp: <2014-07-08 16:40:56 noaki>
 */
#ifndef DIRAC_WILSON_INCLUDED
#define DIRAC_WILSON_INCLUDED

#include "wilsonLikeUtils.hpp"
#include "dirac_WilsonLike.hpp"
#include "Geometry/siteIndex_EvenOdd.hpp"
#include "Geometry/siteIndex.hpp"
#include "Geometry/siteMap.hpp"
#include "eoUtils.hpp"

#ifdef IBM_BGQ_WILSON
#include "bgqwilson.h"
#endif

class Dirac_Wilson: public DiracWilsonLike{
private:
  int Nvol_;
  int Nx_,Ny_,Nz_,Nt_;
  double kpp_;
  GammaMatrix dm_;
  //int boundary[Ndim]; // this is temporary setting.
  const Communicator* comm_;
  const Field* const u_;
  const ffmt_t ff_;
  const gfmt_t gf_;
  
  const std::vector<int> gp_;
  const std::vector<int> gm_;

  EOtype eotype_;
  void init_mult_pm();

  void mult_xp(Field&,const Field&)const;
  void mult_yp(Field&,const Field&)const;
  void mult_zp(Field&,const Field&)const;
  void mult_tp(Field&,const Field&)const;

  void mult_xm(Field&,const Field&)const;
  void mult_ym(Field&,const Field&)const;
  void mult_zm(Field&,const Field&)const;
  void mult_tm(Field&,const Field&)const;

#ifdef IBM_BGQ_WILSON
  void multEO_xp(Field&,const Field&)const;
  void multEO_yp(Field&,const Field&)const;
  void multEO_zp(Field&,const Field&)const;
  void multEO_tp(Field&,const Field&)const;

  void multEO_xm(Field&,const Field&)const;
  void multEO_ym(Field&,const Field&)const;
  void multEO_zm(Field&,const Field&)const;
  void multEO_tm(Field&,const Field&)const;
#endif

  void mult_full(Field&,const Field&)const;   /*! @brief  (1-kpp*D)*f */
  void mult_offdiag(Field&,const Field&)const;/*! @brief  -kpp*D*f */

  // full site operation
  int xsl(int x,int n,int dir)const{return SiteMap::shiftSite.xslice(x,n,dir);}
  int slsize(int x,int dir)const{return SiteMap::shiftSite.slice_size(x,dir);}
  // e/o site operation
  int xsl_e(int x,int n,int dir)const{
    return SiteMap::shiftSite_eo.xslice(x,n,dir);}
  int xsl_o(int x,int n,int dir)const{
    return SiteMap::shiftSite_oe.xslice(x,n,dir);}
  int slsize_e(int x,int dir)const{
    return SiteMap::shiftSite_eo.slice_size(x,dir);}
  int slsize_o(int x,int dir)const{
    return SiteMap::shiftSite_oe.slice_size(x,dir);}

  void mkfrc(Field& fce,const Field& eta,const Field& zeta);

  int(Dirac_Wilson::*slice_out)(int,int,int)const;
  int(Dirac_Wilson::*slice_in)(int,int,int)const;
  int(Dirac_Wilson::*slice_osize)(int,int)const;
  int(Dirac_Wilson::*slice_isize)(int,int)const;

  void(Dirac_Wilson::*mult_core)(Field&,const Field&)const;
  Dirac_Wilson(const Dirac_Wilson&); /*!< @brief simple copy is prohibited.*/

  typedef void(Dirac_Wilson::*mult_pm)(Field&,const Field&)const;

public:
  /*! @brief constructor for e/o-site indexing (odd->even) */
  Dirac_Wilson(double mass,const Field* u,Dop::EOtag)
    :Nvol_(CommonPrms::Nvol()/2),
     Nx_(CommonPrms::Nx()),Ny_(CommonPrms::Ny()),
     Nz_(CommonPrms::Nz()),Nt_(CommonPrms::Nt()),
     kpp_(0.5/(4.0+mass)),ff_(Nvol_),gf_(Nvol_*2),u_(u),
     gp_(SiteIndex_EvenOdd::instance()->esec()),
     gm_(SiteIndex_EvenOdd::instance()->osec()),
     comm_(Communicator::instance()),
     slice_out(&Dirac_Wilson::xsl_e),slice_osize(&Dirac_Wilson::slsize_e),
     slice_in(&Dirac_Wilson::xsl_o),slice_isize(&Dirac_Wilson::slsize_o),
     mult_core(&Dirac_Wilson::mult_offdiag),eotype_(Deo)
  { init_mult_pm(); }

  /*! @brief constructor for e/o-site indexing (even->odd)*/
  Dirac_Wilson(double mass,const Field* u,Dop::OEtag)
    :Nvol_(CommonPrms::Nvol()/2),
     Nx_(CommonPrms::Nx()),Ny_(CommonPrms::Ny()),
     Nz_(CommonPrms::Nz()),Nt_(CommonPrms::Nt()),
     kpp_(0.5/(4.0+mass)),ff_(Nvol_),gf_(Nvol_*2),u_(u),
     gp_(SiteIndex_EvenOdd::instance()->osec()),
     gm_(SiteIndex_EvenOdd::instance()->esec()),
     comm_(Communicator::instance()),
     slice_out(&Dirac_Wilson::xsl_o),slice_osize(&Dirac_Wilson::slsize_o),
     slice_in(&Dirac_Wilson::xsl_e),slice_isize(&Dirac_Wilson::slsize_e),
     mult_core(&Dirac_Wilson::mult_offdiag),eotype_(Doe)
  { init_mult_pm(); }

  ///// manual constructor (regular site-indexing) 
  Dirac_Wilson(double mass,const Field* u)
    :Nvol_(CommonPrms::Nvol()),
     Nx_(CommonPrms::Nx()),Ny_(CommonPrms::Ny()),
     Nz_(CommonPrms::Nz()),Nt_(CommonPrms::Nt()),
     kpp_(0.5/(4.0+mass)),ff_(Nvol_),gf_(Nvol_),u_(u),
     gp_(SiteIndex::instance()->get_lsite()),
     gm_(SiteIndex::instance()->get_lsite()),
     comm_(Communicator::instance()),
     slice_out(&Dirac_Wilson::xsl),slice_osize(&Dirac_Wilson::slsize),
     slice_in(&Dirac_Wilson::xsl),slice_isize(&Dirac_Wilson::slsize),
     mult_core(&Dirac_Wilson::mult_full),eotype_(Void)
  { init_mult_pm(); }

  ///// XML-driven constructor (regular site-indexing) 
  Dirac_Wilson(const XML::node& node,const Field* u)
    :Nvol_(CommonPrms::Nvol()),
     Nx_(CommonPrms::Nx()),Ny_(CommonPrms::Ny()),
     Nz_(CommonPrms::Nz()),Nt_(CommonPrms::Nt()),
     ff_(Nvol_),gf_(Nvol_),u_(u),
     gp_(SiteIndex::instance()->get_lsite()),
     gm_(SiteIndex::instance()->get_lsite()),
     comm_(Communicator::instance()),
     slice_out(&Dirac_Wilson::xsl),slice_osize(&Dirac_Wilson::slsize),
     slice_in(&Dirac_Wilson::xsl),slice_isize(&Dirac_Wilson::slsize),
     mult_core(&Dirac_Wilson::mult_full),eotype_(Void)
  {
    init_mult_pm();
    double mass;
    XML::read(node, "mass", mass,MANDATORY);
    kpp_= 0.5/(4.0+mass);
  }

  mult_pm mult_p[NDIM_];
  mult_pm mult_m[NDIM_];
  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;

#ifdef IBM_BGQ_WILSON
  void mult_ptr(double*,double* const )const;
  void mult_dag_ptr(double*,double* const )const;  
  void mult_ptr_EO(double*,double* const )const;
  void mult_dag_ptr_EO(double*,double* const )const;  
#endif

  const Field gamma5(const Field&) const;
  const Field md_force(const Field&,const Field&)const;
  void md_force_p(Field&,const Field&,const Field&)const;
  void md_force_m(Field&,const Field&,const Field&)const;
  void mkfrc(Field& fce,const Field& eta,const Field& zeta,int mu)const;
  double getKappa() const {return kpp_;}  
  double getMass() const {return 0.5/kpp_-4.0;}  

  size_t fsize()const{return ff_.size();}
  size_t gsize()const{return gf_.size();}
  const Field* getGaugeField_ptr()const{ return u_; }

  void update_internal_state(){}
};
#endif
