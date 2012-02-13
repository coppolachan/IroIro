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

//#define IMPROVED_WILSON

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

private:
  const Field* const u_;
  double kpp_;
  int Nvol_;
  int Ndim_;
  int Nd_;
  //int boundary[Ndim];
  // this is temporary setting.

  mutable ffmt_t* ff_;
  mutable gfmt_t* gf_;

  std::vector<ShiftField*> sf_up_; 
  std::vector<ShiftField*> sf_dn_; 

  const size_t fsize_;
  const size_t gsize_;

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

#ifdef IMPROVED_WILSON
  int r0(int c)const{return 2*c;}
  int r1(int c)const{return 2*(NC_+c);}
  int r2(int c)const{return 2*(2*NC_+c);}
  int r3(int c)const{return 2*(3*NC_+c);} 

  int i0(int c)const{return 2*c+1;}
  int i1(int c)const{return 2*(NC_+c)+1;}
  int i2(int c)const{return 2*(2*NC_+c)+1;}
  int i3(int c)const{return 2*(3*NC_+c)+1;} 
#endif

  int r(int c1,int c2)const{return 2*(NC_*c1+c2);}
  int i(int c1,int c2)const{return 2*(NC_*c1+c2)+1;}

  void mult_full(Field&,const Field&)const;   /*! @brief  -kpp*D*f */
  void mult_offdiag(Field&,const Field&)const;/*! @brief  (1-kpp*D)*f */

  int gsite(int site)const {return site;}
  int esec(int site)const {return SiteIndex_eo::instance()->esec(site);}
  int osec(int site)const {return SiteIndex_eo::instance()->osec(site);}

  int(Dirac_Wilson::*gp)(int)const;
  int(Dirac_Wilson::*gm)(int)const;
  void(Dirac_Wilson::*mult_core)(Field&,const Field&)const;

  Dirac_Wilson(const Dirac_Wilson&); /*!< @brief simple copy is prohibited.*/

public:
  /*! @brief constructor to create instance with e/o site indexing */
  Dirac_Wilson(double mass,const Field* u,Dw::EOtag)
    :kpp_(0.5/(4.0+mass)),u_(u),
     Nvol_(CommonPrms::instance()->Nvol()/2),
     Ndim_(CommonPrms::instance()->Ndim()),
     Nd_(CommonPrms::instance()->Nd()),
     gp(&Dirac_Wilson::esec),
     gm(&Dirac_Wilson::osec),
     mult_core(&Dirac_Wilson::mult_offdiag),
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
     Nd_(CommonPrms::instance()->Nd()),
     gp(&Dirac_Wilson::osec),
     gm(&Dirac_Wilson::esec),
     mult_core(&Dirac_Wilson::mult_offdiag),
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
  /*! @brief constructor to create instance with normal site indexing
   */
  Dirac_Wilson(double mass,const Field* u)
    :kpp_(0.5/(4.0+mass)),u_(u),
     Nvol_(CommonPrms::instance()->Nvol()),
     Ndim_(CommonPrms::instance()->Ndim()),
     Nd_(CommonPrms::instance()->Nd()),
     gp(&Dirac_Wilson::gsite),
     gm(&Dirac_Wilson::gsite),
     mult_core(&Dirac_Wilson::mult_full),
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
     Nd_(CommonPrms::instance()->Nd()),
     gp(&Dirac_Wilson::gsite),
     gm(&Dirac_Wilson::gsite),
     mult_core(&Dirac_Wilson::mult_full),
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
  const ffmt_t get_fermionFormat() const {return *ff_;}
  const std::vector<int> get_gsite() const;

  void update_internal_state(){}
};
#endif
