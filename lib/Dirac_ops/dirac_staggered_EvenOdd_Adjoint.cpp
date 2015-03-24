//--------------------------------------------------------------
/*!@file dirac_wilson_EvenOdd_Adjoint.cpp
  @brief Definition of the Even Odd staggered operator for adjoint rep.
*/
//-------------------------------------------------------------

//Compile only for NC=3
#include "include/macros.hpp"
#if (NC_==3)

#include "dirac_staggered_EvenOdd_Adjoint.hpp"
#include "Fields/field_expressions.hpp"
#include "Tools/fieldUtils.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/sunAdjVec.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include <iostream>

using namespace std;
using namespace FieldUtils;
using namespace SUNmatUtils;
using namespace SUNvecUtils;
using namespace Mapping;

void (Dirac_staggered_EvenOdd_Adjoint::*Dirac_staggered_EvenOdd_Adjoint::
      mult_type[])(Field&,const Field&) const 
= {&Dirac_staggered_EvenOdd_Adjoint::mult_DdagDee,
   &Dirac_staggered_EvenOdd_Adjoint::mult_DdagDoo,
   &Dirac_staggered_EvenOdd_Adjoint::mult_Dfull,};

void (Dirac_staggered_EvenOdd_Adjoint::*Dirac_staggered_EvenOdd_Adjoint::
      mult_dag_type[])(Field&,const Field&) const 
= {&Dirac_staggered_EvenOdd_Adjoint::mult_DdagDee,
   &Dirac_staggered_EvenOdd_Adjoint::mult_DdagDoo,
   &Dirac_staggered_EvenOdd_Adjoint::mult_Dfull_dag,};

void Dirac_staggered_EvenOdd_Adjoint::set_ustag(){
  for(int mu=0; mu<Ndim_; ++mu){
    for(int hs=0; hs<Nvh_; ++hs){
      ue_.data.set(gf_.islice(hs,mu),
		   kse_[Nvh_*mu+hs]
		   *adjoint(SUNmat((*u_)[gff_.islice(idx_->esec(hs),mu)]))); 
      uo_.data.set(gf_.islice(hs,mu),
		   kso_[Nvh_*mu+hs]
		   *adjoint(SUNmat((*u_)[gff_.islice(idx_->osec(hs),mu)]))); 
    }
  }
  ue_*= 0.5/mq_;
  uo_*= 0.5/mq_;
  if(bdry_) bdry_->apply_bc(ue_,uo_);
}

#ifdef IBM_BGQ_WILSON 
#if __VECTOR4DOUBLE__
#include "Architecture_Optimized/dirac_staggered_EvenOdd_Adjoint_BGQ.code"
//DEBUGGING
//#include "dirac_staggered_EvenOdd_Adjoint_improved.code" 
#else
#include "dirac_staggered_EvenOdd_Adjoint_improved.code"
#endif 
#else
#include "dirac_staggered_EvenOdd_Adjoint_improved.code"
#endif 

const Field Dirac_staggered_EvenOdd_Adjoint::mult_eo_dag(const Field& fe) const{
  return -Field(mult_oe(fe));
}

const Field Dirac_staggered_EvenOdd_Adjoint::mult_oe_dag(const Field& fo) const{
  return -Field(mult_eo(fo));
}

const Field Dirac_staggered_EvenOdd_Adjoint::mult(const Field& f)const{
  Field w(fsize_);
  (this->*mult_core)(w,f);
  return w;
}

const Field Dirac_staggered_EvenOdd_Adjoint::mult_dag(const Field& f)const{
  Field w(fsize_);
  (this->*mult_dag_core)(w,f);
  return w;
}

void Dirac_staggered_EvenOdd_Adjoint::mult_DdagDee(Field& we,const Field& fe)const{
  using namespace FieldExpression;
  we = fe -mult_eo(mult_oe(fe));
}

void Dirac_staggered_EvenOdd_Adjoint::mult_DdagDoo(Field& wo,const Field& fo)const{
  using namespace FieldExpression;
  wo = fo -mult_oe(mult_eo(fo));
}

void Dirac_staggered_EvenOdd_Adjoint::mult_Dfull(Field& w,const Field& f)const{
  assert(f.size()==2*fsize_);
  assert(w.size()==2*fsize_);
  w = f;
  w.add(sle_,mult_eo(Field(f[slo_])).getva());
  w.add(slo_,mult_oe(Field(f[sle_])).getva());
}

void Dirac_staggered_EvenOdd_Adjoint::mult_Dfull_dag(Field& w,const Field& f) const{
  assert(f.size()==2*fsize_);
  assert(w.size()==2*fsize_);
  w = f;
  w.add(sle_,-(mult_eo(Field(f[slo_]))).getva());
  w.add(slo_,-(mult_oe(Field(f[sle_]))).getva());
}

const Field Dirac_staggered_EvenOdd_Adjoint::
md_force(const Field& eta,const Field& zeta) const{
  Field fce(gsize_);


  for(int mu=0; mu<Ndim_; ++mu){
    Field etah(fsize_),zetah(fsize_);
    multPoe( etah, eta,mu);
    multPeo(zetah,zeta,mu);

    for(int a=0;a<Nin_;++a){
      for(int b=0;b<Nin_;++b){
	SUNmat Lmd_ab = lmd_commutator(a,b);

	for(int hs=0; hs<Nvh_; ++hs){	
	  fce.add(gff_.islice(idx_->esec(hs),mu),
		  -(Lmd_ab*eta[ff_.index(a,hs)]*zetah[ff_.index(b,hs)]).getva());
	  fce.add(gff_.islice(idx_->osec(hs),mu),
		  (Lmd_ab*zeta[ff_.index(a,hs)]*etah[ff_.index(b,hs)]).getva());
	} 
      }
    }
  }
  fce *= -0.5;
  return fce;
}


#endif //NC=3
