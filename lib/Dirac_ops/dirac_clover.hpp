/*! 
  @file dirac_clover.hpp

  @brief Dirac clover operator class. Class declaration
*/

#ifndef DIRAC_CLOVER_INCLUDED
#define DIRAC_CLOVER_INCLUDED

#include "include/format_F.h"
#include "include/format_G.h"
#include "Dirac_ops/dirac_wilson.hpp"
#include "include/common_fields.hpp"

typedef Format::Format_F ffmt_t;
typedef Format::Format_G gfmt_t;

/*! 
  @ brief Class for the Clover %Dirac operator 
*/
class Dirac_Clover: public DiracWilsonLike{
 
private:
  double csw_;
  const int Nvol_;

  const Field* const u_;
  const Dirac_Wilson* Dw;
  mutable gfmt_t gf_;
  mutable ffmt_t ff_;

  const size_t fsize_;
  const size_t gsize_;

  int gsite(int site)const{return site;}

  static void (Dirac_Clover::*isigma[])(FermionField&,const FermionField&)const;

  void mult_sw(FermionField&, const FermionField&)const ;
  void set_csw() ;
  void set_fieldstrength(GaugeField1D&, const int, const int);

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
  Dirac_Clover(double mass,double csw,const Field* u)
    :csw_(csw),
     Nvol_(CommonPrms::instance()->Nvol()),
     u_(u),
     Dw(new Dirac_Wilson(mass,u)),
     ff_(Nvol_),fsize_(ff_.size()),
     gf_(Nvol_),gsize_(gf_.size()){
    set_csw();
  }

  Dirac_Clover(const XML::node& node,const Field* u)
    :Nvol_(CommonPrms::instance()->Nvol()),
     u_(u),
     Dw(new Dirac_Wilson(node,u)),
     ff_(Nvol_),fsize_(ff_.size()),
     gf_(Nvol_),gsize_(gf_.size()){
    XML::read(node, "Csw", csw_,MANDATORY);
    set_csw();
  }
  
  virtual ~Dirac_Clover(){
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
