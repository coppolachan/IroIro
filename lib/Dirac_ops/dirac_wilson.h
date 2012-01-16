//----------------------------------------------------------------------
// dirac_wilson.h
//----------------------------------------------------------------------
#ifndef DIRAC_WILSON_INCLUDED
#define DIRAC_WILSON_INCLUDED

#include "include/format_F.h"
#include "include/format_G.h"
#include "Main/Geometry/shiftField_eo.h"
#include "Tools/sunMat.h"
#include "Tools/sunVec.h"
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

namespace Dw{
  struct EOtag{};    
  struct OEtag{};    

  double read_mass(const XML::node& node);
}

class Dirac_Wilson: public DiracWilsonLike{
  friend class Dirac_Wilson_EvenOdd;

private:
  const Field* const u_;
  double kpp_;
  int Nvol_;
  int Ndim_;
  //int boundary[Ndim];
  // this is temporary setting.

  mutable ffmt_t* ff_;
  mutable gfmt_t* gf_;

  std::vector<ShiftField*> sf_up_; 
  std::vector<ShiftField*> sf_dn_; 

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

  SUNmat u(int site,int dir) const{
    return SUNmat ((*u_)[gf_->cslice(0,site,dir)]);
  }
  SUNmat u_dag(int site,int dir) const{
    return SUNmat((*u_)[gf_->cslice(0,site,dir)]).dag();
  }

  void mult_xp(Field&,const Field&)const;
  void mult_yp(Field&,const Field&)const;
  void mult_zp(Field&,const Field&)const;
  void mult_tp(Field&,const Field&)const;

  void mult_xm(Field&,const Field&)const;
  void mult_ym(Field&,const Field&)const;
  void mult_zm(Field&,const Field&)const;
  void mult_tm(Field&,const Field&)const;
  void mult_core(Field&,const Field&)const;

  static void(Dirac_Wilson::*mult_p[])(Field&,const Field&)const;
  static void(Dirac_Wilson::*mult_m[])(Field&,const Field&)const;
  int gsite(int site)const {return site;}
  int esec(int site)const {return SiteIndex_eo::instance()->esec(site);}
  int osec(int site)const {return SiteIndex_eo::instance()->osec(site);}

  int(Dirac_Wilson::*gp)(int)const;
  int(Dirac_Wilson::*gm)(int)const;

  void md_force_p(Field&,const Field& eta,const Field& zeta)const;
  void md_force_m(Field&,const Field& eta,const Field& zeta)const;

  /*! @brief private constructor to create instance with e/o site indexing */
  Dirac_Wilson(double mass,const Field* u,Dw::EOtag)
    :kpp_(0.5/(4.0+mass)),u_(u),
     Nvol_(CommonPrms::instance()->Nvol()/2),
     Ndim_(CommonPrms::instance()->Ndim()),
     gp(&Dirac_Wilson::esec),
     gm(&Dirac_Wilson::osec),
     ff_(new ffmt_t(Nvol_)),
     gf_(new gfmt_t(2*Nvol_)),
     fsize_(ff_->size()),
     gsize_(gf_->size()){
    //
    for(int d=0;d<Ndim_;++d){
      sf_up_.push_back(new shift_oup(ff_,d));
      sf_dn_.push_back(new shift_odn(ff_,d));
    }
    CCIO::cout<<"Dirac_Wilson EO created"<<std::endl;
  }

  Dirac_Wilson(double mass,const Field* u,Dw::OEtag)
    :kpp_(0.5/(4.0+mass)),u_(u),
     Nvol_(CommonPrms::instance()->Nvol()/2),
     Ndim_(CommonPrms::instance()->Ndim()),
     gp(&Dirac_Wilson::osec),
     gm(&Dirac_Wilson::esec),
     ff_(new ffmt_t(Nvol_)),
     gf_(new gfmt_t(2*Nvol_)),
     fsize_(ff_->size()),
     gsize_(gf_->size()){
    //
    for(int d=0;d<Ndim_;++d){
      sf_up_.push_back(new shift_eup(ff_,d));
      sf_dn_.push_back(new shift_edn(ff_,d));
    }
    CCIO::cout<<"Dirac_Wilson EO created"<<std::endl;
  }

public:
  /*! @brief constructor to create instance with normal site indexing
   */
  Dirac_Wilson(double mass,const Field* u)
    :kpp_(0.5/(4.0+mass)),u_(u),
     Nvol_(CommonPrms::instance()->Nvol()),
     Ndim_(CommonPrms::instance()->Ndim()),
     gp(&Dirac_Wilson::gsite),
     gm(&Dirac_Wilson::gsite),
     ff_(new ffmt_t(Nvol_)),
     gf_(new gfmt_t(Nvol_)),
     fsize_(ff_->size()),
     gsize_(gf_->size()){
    //
    for(int d=0;d<Ndim_;++d){
      sf_up_.push_back(new shift_up(ff_,d));
      sf_dn_.push_back(new shift_dn(ff_,d));
    }
  }

  Dirac_Wilson(const XML::node& node,const Field* u)
    :u_(u),
     Nvol_(CommonPrms::instance()->Nvol()),
     Ndim_(CommonPrms::instance()->Ndim()),
     gp(&Dirac_Wilson::gsite),
     gm(&Dirac_Wilson::gsite),
     ff_(new ffmt_t(Nvol_)),
     gf_(new gfmt_t(Nvol_)),
     fsize_(ff_->size()),
     gsize_(gf_->size()){
    //
    double mass;
    XML::read(node, "mass", mass);
    kpp_= 0.5/(4.0+mass);
    
    for(int d=0;d<Ndim_;++d){
      sf_up_.push_back(new shift_up(ff_,d));
      sf_dn_.push_back(new shift_dn(ff_,d));
    }
  }

  virtual ~Dirac_Wilson(){
    delete ff_;
    delete gf_;
    for(int d=0;d<sf_up_.size();++d) delete sf_up_[d];
    for(int d=0;d<sf_dn_.size();++d) delete sf_dn_[d];
  }
  
  size_t fsize() const{return fsize_;}
  size_t gsize() const{return gsize_;}

  const Field operator()(int OpType, const Field&)const{};

  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;

  //Preconditioned versions
  const Field mult_prec(const Field& f)const {return f;}//empty now
  const Field mult_dag_prec(const Field& f)const{return f;}//empty now
  const Field left_precond(const Field& f)const{return f;}//empty now
  const Field right_precond(const Field& f)const{return f;}//empty now

  const Field gamma5(const Field&) const;
  const Field proj_p(const Field&) const;
  const Field proj_m(const Field&) const;

  const Field md_force(const Field& eta,const Field& zeta)const;
  
  const double getKappa() const {return kpp_;}  
  const ffmt_t get_fermionFormat() const {return *ff_;}
  const std::vector<int> get_gsite() const;
};
#endif
