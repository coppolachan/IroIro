/*! @file dirac_wilson.hpp
 * @brief Dirac_Wilson class 
 Time-stamp: <2013-04-26 19:47:03 noaki>
 */
#ifndef DIRAC_WILSON_INCLUDED
#define DIRAC_WILSON_INCLUDED

#include "dirac_WilsonLike.hpp"
#include "Geometry/siteIndex_EvenOdd.hpp"
#include "Geometry/siteIndex.hpp"
#include "Geometry/siteMap.hpp"

#ifdef IBM_BGQ_WILSON
#include "bgqwilson.h"
#endif

class Dirac_Wilson: public DiracWilsonLike{

private:
  int Nvol_;
  int Nx_,Ny_,Nz_,Nt_;
  double kpp_;

  //int boundary[Ndim]; // this is temporary setting.
  const Communicator* comm_;
  const Field* const u_;
  const ffmt_t ff_;
  const gfmt_t gf_;

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

  void mult_full(Field&,const Field&)const;   /*! @brief  (1-kpp*D)*f */
  void mult_offdiag(Field&,const Field&)const;/*! @brief  -kpp*D*f */

  int gsite(int site)const {return site;}
  int esec(int hs)const {return SiteIndex_EvenOdd::instance()->esec(hs);}
  int osec(int hs)const {return SiteIndex_EvenOdd::instance()->osec(hs);}

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

  int(Dirac_Wilson::*slice_out)(int,int,int)const;
  int(Dirac_Wilson::*slice_in)(int,int,int)const;
  int(Dirac_Wilson::*slice_osize)(int,int)const;
  int(Dirac_Wilson::*slice_isize)(int,int)const;

  int(Dirac_Wilson::*gp)(int)const;
  int(Dirac_Wilson::*gm)(int)const;

  void(Dirac_Wilson::*mult_core)(Field&,const Field&)const;

  int EO_BGWilson;

  Dirac_Wilson(const Dirac_Wilson&); /*!< @brief simple copy is prohibited.*/

public:
  /*! @brief constructor to create instance with e/o site indexing */
  Dirac_Wilson(double mass,const Field* u,Dop::EOtag)
    :Nvol_(CommonPrms::instance()->Nvol()/2),
     Nx_(CommonPrms::instance()->Nx()),
     Ny_(CommonPrms::instance()->Ny()),
     Nz_(CommonPrms::instance()->Nz()),
     Nt_(CommonPrms::instance()->Nt()),
     kpp_(0.5/(4.0+mass)),ff_(Nvol_),gf_(Nvol_*2),u_(u),
     gp(&Dirac_Wilson::esec),gm(&Dirac_Wilson::osec),
     comm_(Communicator::instance()),
     slice_out(&Dirac_Wilson::xsl_e),slice_osize(&Dirac_Wilson::slsize_e),
     slice_in(&Dirac_Wilson::xsl_o),slice_isize(&Dirac_Wilson::slsize_o),
     mult_core(&Dirac_Wilson::mult_offdiag),
     EO_BGWilson(1){}

  Dirac_Wilson(double mass,const Field* u,Dop::OEtag)
    :Nvol_(CommonPrms::instance()->Nvol()/2),
     Nx_(CommonPrms::instance()->Nx()),
     Ny_(CommonPrms::instance()->Ny()),
     Nz_(CommonPrms::instance()->Nz()),
     Nt_(CommonPrms::instance()->Nt()),
     kpp_(0.5/(4.0+mass)),ff_(Nvol_),gf_(Nvol_*2),u_(u),
     gp(&Dirac_Wilson::osec),gm(&Dirac_Wilson::esec),
     comm_(Communicator::instance()),
     slice_out(&Dirac_Wilson::xsl_o),slice_osize(&Dirac_Wilson::slsize_o),
     slice_in(&Dirac_Wilson::xsl_e),slice_isize(&Dirac_Wilson::slsize_e),
     mult_core(&Dirac_Wilson::mult_offdiag),
     EO_BGWilson(2){}

  /*! @brief constructor to create instance with normal site indexing */
  Dirac_Wilson(double mass,const Field* u)
    :Nvol_(CommonPrms::instance()->Nvol()),
     Nx_(CommonPrms::instance()->Nx()),
     Ny_(CommonPrms::instance()->Ny()),
     Nz_(CommonPrms::instance()->Nz()),
     Nt_(CommonPrms::instance()->Nt()),
     kpp_(0.5/(4.0+mass)),ff_(Nvol_),gf_(Nvol_),u_(u),
     gp(&Dirac_Wilson::gsite),gm(&Dirac_Wilson::gsite),
     comm_(Communicator::instance()),
     slice_out(&Dirac_Wilson::xsl),slice_osize(&Dirac_Wilson::slsize),
     slice_in(&Dirac_Wilson::xsl),slice_isize(&Dirac_Wilson::slsize),
     mult_core(&Dirac_Wilson::mult_full){}

  Dirac_Wilson(const XML::node& node,const Field* u)
    :Nvol_(CommonPrms::instance()->Nvol()),
     Nx_(CommonPrms::instance()->Nx()),
     Ny_(CommonPrms::instance()->Ny()),
     Nz_(CommonPrms::instance()->Nz()),
     Nt_(CommonPrms::instance()->Nt()),
     ff_(Nvol_),gf_(Nvol_),u_(u),
     gp(&Dirac_Wilson::gsite),gm(&Dirac_Wilson::gsite),
     comm_(Communicator::instance()),
     slice_out(&Dirac_Wilson::xsl),slice_osize(&Dirac_Wilson::slsize),
     slice_in(&Dirac_Wilson::xsl),slice_isize(&Dirac_Wilson::slsize),
     mult_core(&Dirac_Wilson::mult_full){
    //
    double mass;
    XML::read(node, "mass", mass);
    kpp_= 0.5/(4.0+mass);
  }

  virtual ~Dirac_Wilson(){}
  
  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;

#ifdef IBM_BGQ_WILSON
  void mult_ptr(double*,double* const )const;
  void mult_dag_ptr(double*,double* const )const;  
  void mult_ptr_EO(double*,double* const )const;
  void mult_dag_ptr_EO(double*,double* const )const;  
#endif

  ////////////////////////////////////////Preconditioned versions
  // Wilson operator has no defined preconditioner now 
  const Field mult_prec     (const Field& f)const{return f;}
  const Field mult_dag_prec (const Field& f)const{return f;}
  const Field left_prec     (const Field& f)const{return f;}
  const Field right_prec    (const Field& f)const{return f;}
  const Field left_dag_prec (const Field& f)const{return f;}
  const Field right_dag_prec(const Field& f)const{return f;}
  //////////////////////////////////////////////////////////////
  const Field gamma5(const Field&) const;

  const Field md_force(const Field&,const Field&)const;
  void md_force_p(Field&,const Field&,const Field&)const;
  void md_force_m(Field&,const Field&,const Field&)const;
  double getKappa() const {return kpp_;}  

  size_t fsize()const{return ff_.size();}
  size_t gsize()const{return gf_.size();}

  void update_internal_state(){}
};
#endif
