/*! 
  @file dirac_clover.hpp

  @brief Dirac clover operator class. Class declaration

*/

#ifndef DIRAC_CLOVER_INCLUDED
#define DIRAC_CLOVER_INCLUDED

#include "dirac.hpp"
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
  mutable gfmt_t* gf_;
  mutable ffmt_t* ff_;

  const size_t fsize_;
  const size_t gsize_;

  int gsite(int site)const {return site;}

  //void (Dirac_Clover::*mult_isigma[])(Field&,const Field&) const;

  void mult_isigma(FermionField&, const FermionField&,int mu,int nu) const;
  
  void mult_csw(FermionField&, const FermionField&) const ;
  void set_csw() ;
  void set_fieldstrength(GaugeField1D&, const int mu, const int nu);

  void mult_isigma23(FermionField&, const FermionField&)const;
  void mult_isigma12(FermionField&, const FermionField&)const;
  void mult_isigma31(FermionField&, const FermionField&)const;
  void mult_isigma41(FermionField&, const FermionField&)const;
  void mult_isigma42(FermionField&, const FermionField&)const;
  void mult_isigma43(FermionField&, const FermionField&)const;

  const Field md_force_block(const FermionField&,const FermionField&)const;

  GaugeField1D d_Bx, d_By, d_Bz, d_Ex, d_Ey, d_Ez;
  // Bx = -iF(1,2), By = -iF(2,1), -iBz = F(0,1)
  // Ex = -iF(4,0), Ey = -iF(4,1), Ez = -iF(4,2)

  //auxiliary, eventually moved outside
  //const std::valarray<double> anti_herm(const SUNmat& m);
  void external_prod(GaugeField1D& res, const FermionField& A, const FermionField& B) const;

public:
  Dirac_Clover(double mass,double csw,const Field* u)
    :csw_(csw),
     Nvol_(CommonPrms::instance()->Nvol()),
     u_(u),
     Dw(new Dirac_Wilson(mass,u)),
     gf_(new gfmt_t(Nvol_)),
     ff_(new ffmt_t(Nvol_)),
     fsize_(ff_->size()),
     gsize_(gf_->size()){
    set_csw();
  }

  Dirac_Clover(const XML::node& node,const Field* u)
    :Nvol_(CommonPrms::instance()->Nvol()),
     u_(u),
     Dw(new Dirac_Wilson(node,u)),
     gf_(new gfmt_t(Nvol_)),
     ff_(new ffmt_t(Nvol_)),
     fsize_(ff_->size()),
     gsize_(gf_->size()){

    XML::read(node, "Csw", csw_,MANDATORY);
    
    set_csw();
  }

  virtual ~Dirac_Clover(){
    delete ff_;
    delete gf_;
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
  const ffmt_t get_fermionFormat() const {return *ff_;}
};
#endif
