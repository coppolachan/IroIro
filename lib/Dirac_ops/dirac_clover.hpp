/*! 
  @file dirac_clover.hpp
  @brief Dirac clover operator class. Class declaration
  Time-stamp:<>
*/

#ifndef DIRAC_CLOVER_INCLUDED
#define DIRAC_CLOVER_INCLUDED

#include "include/format_F.h"
#include "include/format_G.h"
#include "include/common_fields.hpp"
#include "dirac_WilsonLike.hpp"
/*! 
  @ brief Class for the Clover Dirac operator 
*/
class Dirac_Clover: public DiracWilsonLike{
 
private:
  const DiracWilsonLike* Dw_;
  const Field* const u_;
  double csw_;
  double kpp_;
  int Nvol_;

  const ffmt_t ff_;
  const gfmt_t gf_;

  static void (Dirac_Clover::*isigma[])(FermionField&,const FermionField&)const;

  void mult_sw(FermionField&, const FermionField&)const ;
  void set_csw() ;
  void set_fieldstrength(GaugeField1D&,int,int);

  void isigma_diag(FermionField&, const FermionField&)const;

  void isigma_12(FermionField&, const FermionField&)const;
  void isigma_13(FermionField&, const FermionField&)const;
  void isigma_14(FermionField&, const FermionField&)const;

  void isigma_21(FermionField&, const FermionField&)const;
  void isigma_23(FermionField&, const FermionField&)const;
  void isigma_24(FermionField&, const FermionField&)const;

  void isigma_31(FermionField&, const FermionField&)const;
  void isigma_32(FermionField&, const FermionField&)const;
  void isigma_34(FermionField&, const FermionField&)const;

  void isigma_41(FermionField&, const FermionField&)const;
  void isigma_42(FermionField&, const FermionField&)const;
  void isigma_43(FermionField&, const FermionField&)const;

  const Field md_force_block(const FermionField&,const FermionField&)const;

  GaugeField1D Bx_, By_, Bz_, Ex_, Ey_, Ez_;
  // Bx = -iF(1,2), By = -iF(2,1), -iBz = F(0,1)
  // Ex = -iF(4,0), Ey = -iF(4,1), Ez = -iF(4,2)

public:
  Dirac_Clover(DiracWilsonLike* Dw,double csw,const Field* u)
    :Nvol_(CommonPrms::Nvol()),Dw_(Dw),u_(u),
     csw_(csw),kpp_(0.5/(Dw->getMass()+4.0)),
     ff_(Nvol_),gf_(Nvol_){ set_csw(); }
  
  Dirac_Clover(const XML::node& node,DiracWilsonLike* Dw,const Field* u)
    :Nvol_(CommonPrms::Nvol()),Dw_(Dw),u_(u),
     kpp_(0.5/(Dw->getMass()+4.0)),
     ff_(Nvol_),gf_(Nvol_){
    //
    XML::read(node, "Csw", csw_,MANDATORY);
    set_csw();
  }
  
  virtual ~Dirac_Clover(){ delete Dw_; }
  
  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;

  size_t fsize() const{return ff_.size();}
  size_t gsize() const{return gf_.size();}
  double getMass()const{return Dw_->getMass();}  

  const Field* getGaugeField_ptr()const{return u_;}

  const Field gamma5(const Field&) const;
  const Field md_force(const Field& eta,const Field& zeta)const;
  void update_internal_state() {set_csw();}
};
#endif
