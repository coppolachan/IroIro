//----------------------------------------------------------------------
// dirac_wilson.hpp
//----------------------------------------------------------------------
#ifndef DIRAC_WILSON_INCLUDED
#define DIRAC_WILSON_INCLUDED

#include "dirac.hpp"
#include "include/format_F.h"
#include "include/format_G.h"
#include "include/pugi_interface.h"

#include "Geometry/siteIndex_EvenOdd.hpp"
#include "Geometry/siteIndex.hpp"
#include "Geometry/siteMap.hpp"

#ifdef IBM_BGQ_WILSON
#include "bgqwilson.h"
#endif

typedef Format::Format_F ffmt_t;
typedef Format::Format_G gfmt_t;

class Dirac_Wilson: public DiracWilsonLike{

private:
  const Field* const u_;
  double kpp_;
  int Nx_,Ny_,Nz_,Nt_,Nvol_;

  //int boundary[Ndim]; // this is temporary setting.

  const ffmt_t ff_;
  const gfmt_t gf_;

  const size_t fsize_;
  const size_t gsize_;
  const Communicator* comm_;

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

  int r0(int c)const{return 2*c;}
  int r1(int c)const{return 2*(NC_+c);}
  int r2(int c)const{return 2*(2*NC_+c);}
  int r3(int c)const{return 2*(3*NC_+c);} 

  int i0(int c)const{return 2*c+1;}
  int i1(int c)const{return 2*(NC_+c)+1;}
  int i2(int c)const{return 2*(2*NC_+c)+1;}
  int i3(int c)const{return 2*(3*NC_+c)+1;} 

  int re(int c1,int c2)const{return 2*(NC_*c1+c2);}
  int im(int c1,int c2)const{return 2*(NC_*c1+c2)+1;}

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
    :kpp_(0.5/(4.0+mass)),u_(u),
     Nx_(CommonPrms::instance()->Nx()),
     Ny_(CommonPrms::instance()->Ny()),
     Nz_(CommonPrms::instance()->Nz()),
     Nt_(CommonPrms::instance()->Nt()),
     Nvol_(CommonPrms::instance()->Nvol()/2),
     gp(&Dirac_Wilson::esec),gm(&Dirac_Wilson::osec),
     comm_(Communicator::instance()),
     slice_out(&Dirac_Wilson::xsl_e),slice_osize(&Dirac_Wilson::slsize_e),
     slice_in(&Dirac_Wilson::xsl_o),slice_isize(&Dirac_Wilson::slsize_o),
     mult_core(&Dirac_Wilson::mult_offdiag),
     ff_(Nvol_),  fsize_(ff_.size()),
     gf_(2*Nvol_),gsize_(gf_.size()),
     EO_BGWilson(1){}

  Dirac_Wilson(double mass,const Field* u,Dop::OEtag)
    :kpp_(0.5/(4.0+mass)),u_(u),
     Nx_(CommonPrms::instance()->Nx()),
     Ny_(CommonPrms::instance()->Ny()),
     Nz_(CommonPrms::instance()->Nz()),
     Nt_(CommonPrms::instance()->Nt()),
     Nvol_(CommonPrms::instance()->Nvol()/2),
     gp(&Dirac_Wilson::osec),gm(&Dirac_Wilson::esec),
     comm_(Communicator::instance()),
     slice_out(&Dirac_Wilson::xsl_o),slice_osize(&Dirac_Wilson::slsize_o),
     slice_in(&Dirac_Wilson::xsl_e),slice_isize(&Dirac_Wilson::slsize_e),
     mult_core(&Dirac_Wilson::mult_offdiag),
     ff_(Nvol_),  fsize_(ff_.size()),
     gf_(2*Nvol_),gsize_(gf_.size()),
     EO_BGWilson(2){}

  /*! @brief constructor to create instance with normal site indexing */
  Dirac_Wilson(double mass,const Field* u)
    :kpp_(0.5/(4.0+mass)),u_(u),
     Nx_(CommonPrms::instance()->Nx()),
     Ny_(CommonPrms::instance()->Ny()),
     Nz_(CommonPrms::instance()->Nz()),
     Nt_(CommonPrms::instance()->Nt()),
     Nvol_(CommonPrms::instance()->Nvol()),
     gp(&Dirac_Wilson::gsite),gm(&Dirac_Wilson::gsite),
     comm_(Communicator::instance()),
     slice_out(&Dirac_Wilson::xsl),slice_osize(&Dirac_Wilson::slsize),
     slice_in(&Dirac_Wilson::xsl),slice_isize(&Dirac_Wilson::slsize),
     mult_core(&Dirac_Wilson::mult_full),
     ff_(Nvol_),fsize_(ff_.size()),
     gf_(Nvol_),gsize_(gf_.size()){}

  Dirac_Wilson(const XML::node& node,const Field* u)
    :u_(u),
     Nx_(CommonPrms::instance()->Nx()),
     Ny_(CommonPrms::instance()->Ny()),
     Nz_(CommonPrms::instance()->Nz()),
     Nt_(CommonPrms::instance()->Nt()),
     Nvol_(CommonPrms::instance()->Nvol()),
     gp(&Dirac_Wilson::gsite),gm(&Dirac_Wilson::gsite),
     comm_(Communicator::instance()),
     slice_out(&Dirac_Wilson::xsl),slice_osize(&Dirac_Wilson::slsize),
     slice_in(&Dirac_Wilson::xsl),slice_isize(&Dirac_Wilson::slsize),
     mult_core(&Dirac_Wilson::mult_full),
     ff_(Nvol_),fsize_(ff_.size()),
     gf_(Nvol_),gsize_(gf_.size()){
    //
    double mass;
    XML::read(node, "mass", mass);
    kpp_= 0.5/(4.0+mass);
  }

  virtual ~Dirac_Wilson(){  }
  
  size_t fsize() const{return fsize_;}
  size_t gsize() const{return gsize_;}

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
  void gamma5_mult(Field&, const Field&) const;
  void gamma5_ptr(double*, double* const) const;
  const Field proj_p(const Field&) const;
  const Field proj_m(const Field&) const;
  void proj_p(Field&, const Field&, int) const;
  void proj_m(Field&, const Field&, int) const;

  const Field md_force(const Field& , const Field&)const;
  void md_force_p(Field&,const Field&,const Field&)const;
  void md_force_m(Field&,const Field&,const Field&)const;
  void get_RandGauss(std::valarray<double>& phi,const RandNum& rng)const;
  Format::Format_F getFermionFormat()const{return ff_;}

  void update_internal_state(){}
  double getKappa() const {return kpp_;}  
  int Nvol()const{return Nvol_;}
};
#endif
