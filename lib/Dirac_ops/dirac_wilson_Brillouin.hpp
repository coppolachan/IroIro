/*! @filename dirac_wilson_Brillouin.hpp
 * @brief declaration of Dirac_Wilson_Brillouin class 
 * Time-stamp: <2013-04-19 17:54:09 noaki>
 ----------------------------------------------------------------------*/
#ifndef DIRAC_WILSON_BRILLOUIN_INCLUDED
#define DIRAC_WILSON_BRILLOUIN_INCLUDED

#include "dirac.hpp"
#include "include/format_F.h"
#include "include/format_G.h"
#include "include/pugi_interface.h"

#include "Main/Geometry/siteIndex.hpp"
#include "Main/Geometry/siteMap.hpp"

typedef Format::Format_F ffmt_t;
typedef Format::Format_G gfmt_t;

class Dirac_Wilson_Brillouin: public DiracWilsonLike{
private:
  const Field* const u_;
  double kbr_,m_;
  int Nx_,Ny_,Nz_,Nt_,Nvol_,Nin_;

  const Communicator* comm_;

  int sr(int s,int c)const{return 2*(s*NC_+c);}
  int si(int s,int c)const{return 2*(s*NC_+c)+1;}
  
  // full site operation
  int xsl(int x,int n,int dir)const{return SiteMap::shiftSite.xslice(x,n,dir);}
  int slsize(int dir)const{return SiteMap::shiftSite.slice_size(0,dir);}

  const Field sft_p(const Field& f,int dir)const;
  const Field sft_m(const Field& f,int dir)const;
  const Field lap(const Field& f,int dir,double a)const;
  const Field del(const Field& f,int dir,double a)const;

  const Field gamma_x(const Field&)const;
  const Field gamma_y(const Field&)const;
  const Field gamma_z(const Field&)const;
  const Field gamma_t(const Field&)const;

  static const Field (Dirac_Wilson_Brillouin::*gm[])(const Field&)const;
  
  int sgm(int,int,int)const;
public:
 /*! @brief constructor to create instance with normal site indexing */
  Dirac_Wilson_Brillouin(double mass,const Field* u)
    :kbr_(1.0/((15.0/8.0)+mass)),u_(u),
     Nx_(CommonPrms::instance()->Nx()),
     Ny_(CommonPrms::instance()->Ny()),
     Nz_(CommonPrms::instance()->Nz()),
     Nt_(CommonPrms::instance()->Nt()),
     Nvol_(CommonPrms::instance()->Nvol()),
     comm_(Communicator::instance()),
     ff_(Nvol_),fsize_(ff_.size()),
     gf_(Nvol_),gsize_(gf_.size()),
     m_(mass),
     Nin_(ff_.Nin()){}

  Dirac_Wilson_Brillouin(const XML::node& node,const Field* u)
    :u_(u),
     Nx_(CommonPrms::instance()->Nx()),
     Ny_(CommonPrms::instance()->Ny()),
     Nz_(CommonPrms::instance()->Nz()),
     Nt_(CommonPrms::instance()->Nt()),
     Nvol_(CommonPrms::instance()->Nvol()),
     comm_(Communicator::instance()),
     ff_(Nvol_),fsize_(ff_.size()),
     gf_(Nvol_),gsize_(gf_.size()),
     Nin_(ff_.Nin()){
    //
    double mass;
    XML::read(node,"mass",mass);
    kbr_= 1.0/((15.0/8.0)+mass);
    m_ = mass;
  }

  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;
  const Field mult_H(const Field&)const;
  const Field mult_del(const Field&)const;
  const Field mult_lap(const Field&)const;

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
  const Field md_force(const Field& , const Field&)const{}
  void md_force_p(Field&,const Field&,const Field&)const{}
  void md_force_m(Field&,const Field&,const Field&)const{}

  double getKappa() const {return kbr_;} 
};

#endif
