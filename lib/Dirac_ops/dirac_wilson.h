//----------------------------------------------------------------------
// dirac_wilson.h
//----------------------------------------------------------------------
#ifndef DIRAC_WILSON_INCLUDED
#define DIRAC_WILSON_INCLUDED

#ifndef FORMAT_F_INCLUDED
#include "include/format_F.h"
#endif

#ifndef FORMAT_G_INCLUDED
#include "include/format_G.h"
#endif

#ifndef SHIFTFIELD_EO_INCLUDED
#include "Main/Geometry/shiftField_eo.h"
#endif

#ifndef SUNMAT_INCLUDED
#include "Tools/sunMat.h"
#endif

#ifndef SUNVEC_INCLUDED
#include "Tools/sunVec.h"
#endif

#include "dirac.h"
#include "include/pugi_interface.h"

typedef Format::Format_F ffmt_t;
typedef Format::Format_G gfmt_t;

typedef ShiftField_up<ffmt_t> shift_up;
typedef ShiftField_dn<ffmt_t> shift_dn;

typedef ShiftField_even_up<ffmt_t> shift_eup;
typedef ShiftField_even_dn<ffmt_t> shift_edn;
typedef ShiftField_odd_up<ffmt_t>  shift_oup;
typedef ShiftField_odd_dn<ffmt_t>  shift_odn;

struct Dirac_WilsonParams{
  double kappa;
  
  Dirac_WilsonParams(const XML::node node){
    double mass;
    XML::read(node, "mass", mass);
    kappa = 0.5/(4.0+mass);
  }
  Dirac_WilsonParams(const double mass){
    kappa = 0.5/(4.0+mass);
  } 


};



class Dirac_Wilson : public DiracWilsonLike {
protected:
  const Field* const u_;
  const Dirac_WilsonParams Params;
  int Nvol_;
  int Ndim_;

  mutable ffmt_t* ff_;
  mutable gfmt_t* gf_;
  ShiftField* sf_up_; 
  ShiftField* sf_dn_; 


  const size_t fsize_;
  size_t gsize_;

  SUNvec v(const Field& f,int spin,int site) const{
    return SUNvec(f[ff_->cslice(spin,site)]);
  }
  SUNvec v(const ShiftField* sf,int spin,int site) const{
    return SUNvec(sf->cv(spin,site));
  }

  SUNvec v_Ix(const Field& f,int spin,int site) const{
    return SUNvec(f[ff_->cslice(spin,site)]).xI();
  }
  SUNvec v_Ix(const ShiftField* sf,int spin,int site) const{
    return SUNvec(sf->cv(spin,site)).xI();
  }
  
  virtual SUNmat u(int site,int dir) const{
    return SUNmat((*u_)[gf_->cslice(0,site,dir)]);
  }
  virtual SUNmat u_dag(int site,int dir) const{
    return SUNmat((*u_)[gf_->cslice(0,site,dir)]).dag();
  }

  void mult_xp(Field&,ShiftField*)const;
  void mult_yp(Field&,ShiftField*)const;
  void mult_zp(Field&,ShiftField*)const;
  void mult_tp(Field&,ShiftField*)const;

  void mult_xm(std::valarray<double>&,const Field&)const;
  void mult_ym(std::valarray<double>&,const Field&)const;
  void mult_zm(std::valarray<double>&,const Field&)const;
  void mult_tm(std::valarray<double>&,const Field&)const;

  static void(Dirac_Wilson::*mult_p[])(Field&,ShiftField*)const;
  static void(Dirac_Wilson::*mult_m[])(std::valarray<double>&,
				       const Field&)const;
  
  void mult_core(Field&, const Field&)const;

public:
  Dirac_Wilson(const double mass,const Field* u)
    :Params(Dirac_WilsonParams(mass)),u_(u),
     Nvol_(CommonPrms::instance()->Nvol()),
     Ndim_(CommonPrms::instance()->Ndim()),
     ff_(new ffmt_t(Nvol_)), 
     gf_(new gfmt_t(Nvol_)),
     sf_up_( new shift_up(ff_)),
     sf_dn_( new shift_dn(ff_)),
     fsize_(ff_->size()),
     gsize_(gf_->size()){}  

  Dirac_Wilson(const XML::node node,const Field* u)
    :Params(Dirac_WilsonParams(node)),u_(u),
     Nvol_(CommonPrms::instance()->Nvol()),
     Ndim_(CommonPrms::instance()->Ndim()),
     ff_(new ffmt_t(Nvol_)), 
     gf_(new gfmt_t(Nvol_)),
     sf_up_( new shift_up(ff_)),
     sf_dn_( new shift_dn(ff_)),
     fsize_(ff_->size()),
     gsize_(gf_->size()){}  


  ~Dirac_Wilson(){
    delete ff_;
    delete gf_;
    delete sf_up_;
    delete sf_dn_;
   }
  
  size_t fsize() const{return fsize_;}
  size_t gsize() const{return gsize_;}

  const Field operator()(int OpType, const Field&)const{}; 

  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;
  const Field gamma5(const Field&) const;
  const Field md_force(const Field& eta,const Field& zeta)const;

  const double getKappa() const {return Params.kappa;}

};

#endif
