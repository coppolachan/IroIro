/*! 
  @file dirac_clover.hpp
  @brief Dirac clover operator class. Class declaration
*/
#ifndef DIRAC_CLOVER_INCLUDED
#define DIRAC_CLOVER_INCLUDED

#include "dirac.h"
#include "include/format_F.h"
#include "include/format_G.h"
#include "Tools/sunVec.h"
#include "Main/Geometry/shiftField.h"
#include "Measurements/GaugeM/staples.hpp"
#include "Dirac_ops/dirac_wilson.h"

typedef Format::Format_F ffmt_t;
typedef Format::Format_G gfmt_t;
typedef ShiftField_up<ffmt_t> shift_up;
typedef ShiftField_dn<ffmt_t> shift_dn;

/*! 
  @ brief Class for the Clover %Dirac operator 
*/
class Dirac_Clover: public DiracWilsonLike{
 
private:
  double csw_;
  const int Nvol_;
  const int Ndim_;

  const Field* const u_;
  const Dirac_Wilson* Dw;
  mutable gfmt_t gf_;
  mutable ffmt_t ff_;
  const Staples stpl_;

  std::vector<ShiftField*> sf_up_; 
  std::vector<ShiftField*> sf_dn_; 

  const size_t fsize_;
  const size_t gsize_;

  SUNvec v(const Field& f,int spin,int site) const{
    return SUNvec(f[ff_.cslice(spin,site)]);
  }
  SUNvec v(const ShiftField* sf,int spin,int site) const{
    return SUNvec(sf->cv(spin,site));
  }

  int gsite(int site)const {return site;}

  static void (Dirac_Clover::*isigma[])(Field&,const Field&)const;

  void mult_sw(Field&, const Field&)const ;
  void set_csw() ;
  void set_fieldstrength(GaugeField1D&, const int mu, const int nu);

  void isigma_diag(Field&, const Field&)const;

  void isigma_12(Field&, const Field&)const;
  void isigma_13(Field&, const Field&)const;
  void isigma_14(Field&, const Field&)const;

  void isigma_21(Field&, const Field&)const;
  void isigma_23(Field&, const Field&)const;
  void isigma_24(Field&, const Field&)const;

  void isigma_31(Field&, const Field&)const;
  void isigma_32(Field&, const Field&)const;
  void isigma_34(Field&, const Field&)const;

  void isigma_41(Field&, const Field&)const;
  void isigma_42(Field&, const Field&)const;
  void isigma_43(Field&, const Field&)const;

  const Field md_force_block(const Field&,const Field&)const;

  GaugeField1D Bx_, By_, Bz_, Ex_, Ey_, Ez_;
  // Bx = -iF(1,2), By = -iF(2,1), -iBz = F(0,1)
  // Ex = -iF(4,0), Ey = -iF(4,1), Ez = -iF(4,2)

public:
  Dirac_Clover(double mass,double csw,const Field* u)
    :csw_(csw),
     Nvol_(CommonPrms::instance()->Nvol()),
     Ndim_(CommonPrms::instance()->Ndim()),
     u_(u),
     Dw(new Dirac_Wilson(mass,u)),
     ff_(Nvol_),fsize_(ff_.size()),
     gf_(Nvol_),gsize_(gf_.size()),
     stpl_(gf_){
    for(int d=0;d<Ndim_;++d){
      sf_up_.push_back(new shift_up(&ff_,d));
      sf_dn_.push_back(new shift_dn(&ff_,d));
    }
    set_csw();
  }

  Dirac_Clover(const XML::node& node,const Field* u)
    :Nvol_(CommonPrms::instance()->Nvol()),
     Ndim_(CommonPrms::instance()->Ndim()),
     u_(u),
     Dw(new Dirac_Wilson(node,u)),
     ff_(Nvol_),fsize_(ff_.size()),
     gf_(Nvol_),gsize_(gf_.size()),
     stpl_(gf_){

    XML::read(node, "Csw", csw_,MANDATORY);
    for(int d=0;d<Ndim_;++d){
      sf_up_.push_back(new shift_up(&ff_,d));
      sf_dn_.push_back(new shift_dn(&ff_,d));
    }
    set_csw();
  }

  virtual ~Dirac_Clover(){
    for(int d=0;d<sf_up_.size();++d) delete sf_up_[d];
    for(int d=0;d<sf_dn_.size();++d) delete sf_dn_[d];
    delete Dw;
  }
  
  size_t fsize() const{return fsize_;}
  size_t gsize() const{return gsize_;}

  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;

  ////////////////////////////////////////Preconditioned versions
  // Clover operator has no preconditioner now 
  const Field mult_prec     (const Field&f)const{return f;}
  const Field mult_dag_prec (const Field&f)const{return f;}
  const Field left_prec     (const Field&f)const{return f;}
  const Field right_prec    (const Field&f)const{return f;}
  const Field left_dag_prec (const Field&f)const{return f;}
  const Field right_dag_prec(const Field&f)const{return f;}
  //////////////////////////////////////////////////////////////

  const Field gamma5(const Field&) const;
  const Field md_force(const Field& eta,const Field& zeta)const;

  void update_internal_state() {set_csw();}
  
  const std::vector<int> get_gsite() const;
  const ffmt_t get_fermionFormat() const {return ff_;}
};
#endif
