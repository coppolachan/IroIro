//--------------------------------------------------------------
/*!@file dirac_wilson_EvenOdd_Adjoint.cpp
  @brief Definition of the Even Odd staggered operator for adjoint rep.
*/
//-------------------------------------------------------------
#include "dirac_staggered_EvenOdd_Adjoint.hpp"
#include "Fields/field_expressions.hpp"
#include "Tools/fieldUtils.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/sunAdjVec.hpp"
#include "Tools/randNum_MP.h"
#include "Measurements/GaugeM/staples.hpp"
#include <iostream>

using namespace std;
using namespace FieldUtils;
using namespace SUNmatUtils;
using namespace SUNvecUtils;
using namespace Mapping;

void Dirac_staggered_EvenOdd_Adjoint::
get_RandGauss(valarray<double>& phi,const RandNum& rng)const{
  MPrand::mp_get(phi,rng,SiteIndex_EvenOdd::instance()->get_gsite(),ff_);
}

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
  /*
  Staples stpl;
  double plaq = stpl.plaquette(GaugeField(*u_));
  CCIO::cout<<"Dirac_staggered_EvenOdd_Adjoint::set_ustag   plaq="<<plaq<<"\n";
  */
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

#if 1
#include "dirac_staggered_EvenOdd_Adjoint_improved.code"
#endif 

#if 0
void Dirac_staggered_EvenOdd_Adjoint::
multPeo(Field& we,const Field& fo,int mu)const{
  AdjFermionField1sp ft = shiftField_eo(AdjFermionField1sp(fo),mu,Forward());
  for(int hs=0; hs<Nvh_; ++hs) 
    we.add(ff_.islice(hs),(mat(ue_,hs,mu)*vec(ft,hs)).getva());
}

void Dirac_staggered_EvenOdd_Adjoint::
multPoe(Field& wo,const Field& fe,int mu)const{
  AdjFermionField1sp ft = shiftField_oe(AdjFermionField1sp(fe),mu,Forward());
  for(int hs=0; hs<Nvh_; ++hs)
    wo.add(ff_.islice(hs),(mat(uo_,hs,mu)*vec(ft,hs)).getva());
}

const Field Dirac_staggered_EvenOdd_Adjoint::mult_eo(const Field& fo)const{
  Field we(fsize_);
  for(int mu=0; mu<Ndim_; ++mu){
    multPeo(we,fo,mu);           //forward differenciation
    
    AdjFermionField1sp ft(Nvh_);    
    for(int hs=0; hs<Nvh_; ++hs) //backward differenciation
      ft.data.set(ff_.islice(hs),
		  (mat_dag(uo_,hs,mu)*SUNadjVec(fo[ff_.islice(hs)])).getva());
    we -= shiftField_eo(ft,mu,Backward()).data;
  }
  return we;
}

const Field Dirac_staggered_EvenOdd_Adjoint::mult_oe(const Field& fe)const{
  Field wo(fsize_);
  for(int mu=0; mu<Ndim_; ++mu){
    multPoe(wo,fe,mu);           //forward differenciation
    
    AdjFermionField1sp ft(Nvh_);    
    for(int hs=0; hs<Nvh_; ++hs) //backward differenciation
      ft.data.set(ff_.islice(hs),
		  (mat_dag(ue_,hs,mu)*SUNadjVec(fe[ff_.islice(hs)])).getva());

    wo -= shiftField_oe(ft,mu,Backward()).data;
  }
  return wo;
}
#endif //0

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
	  std::slice xsl = ff_.islice(hs);
	  fce.add(gff_.islice(idx_->esec(hs),mu),
		  -(Lmd_ab*eta[ff_.index(a,hs)]*zetah[ff_.index(b,hs)]).getva());
	  fce.add(gff_.islice(idx_->osec(hs),mu),
		  (Lmd_ab*zeta[ff_.index(a,hs)]*etah[ff_.index(b,hs)]).getva());
	} 
      }
    }
  }
  fce *= 0.5;
  return fce;
}
