/*! @file cloverTerm.hpp
 *  @brief operation of the clover term 
 *  Time-stamp: <2014-08-29 18:44:53 noaki>
 */
#ifndef CLOVERTERM_INCLUDED
#define CLOVERTERM_INCLUDED

#include "include/common_fields.hpp"
#include "wilsonLikeUtils.hpp"
#include "dirac_WilsonLike.hpp"

class CloverTerm {
private:
  const Field* const u_;
  double csw_;
  int Nvol_;
  const ffmt_t ff_;
  const gfmt_t gf_;
  GammaMatrix dm_;

  //static void (CloverTerm::*isigma[])(FermionField&,const FermionField&)const;
  //void isigma_diag(FermionField&, const FermionField&)const;

  void mult_sw(Field&, const Field&)const;

  mutable GaugeField1D Bx_, By_, Bz_, Ex_, Ey_, Ez_;
  // Bx = -iF(1,2), By = -iF(2,1), -iBz = F(0,1)
  // Ex = -iF(4,0), Ey = -iF(4,1), Ez = -iF(4,2)

  void mat_vec(double* wp,double* Fp,double* vp)const;
  int re(int c1,int c2)const{return 2*(NC_*c1+c2);}
  int im(int c1,int c2)const{return 2*(NC_*c1+c2)+1;}
  int fr(int s,int c)const{return 2*(NC_*s+c);}
  int fi(int s,int c)const{return 2*(NC_*s+c)+1;}

public:
  CloverTerm(double csw,const Field* u)
    :Nvol_(CommonPrms::Nvol()),csw_(csw),
     u_(u),ff_(Nvol_),gf_(Nvol_,1){ set_csw();}
  
  const Field mult(const Field&)const;
  void set_csw();
};

#endif
