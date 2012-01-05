//--------------------------------------------------------------
/*
  @file dirac_wilson_EvenOdd.cpp
  @brief Definition of Even Odd wilson operator
*/
//--------------------------------------------------------------
#include "dirac_wilson_EvenOdd.h"
using namespace std;

const Field Dw::EvOd::mult(const Field& f)const{
  Field w(f.size());
  Dirac_Wilson::mult_core(w,const_cast<Field&>(f));
  return w;
}

const Field Dw::EvOd::mult_dag(const Field& f)const{
  return gamma5(mult(gamma5(f)));
}
const Field Dw::EvOd::md_force(const Field& eta,const Field& zeta) const{
  return Dirac_Wilson::md_force(eta,zeta);
}

const vector<int> Dw::EvOd::get_gsite() const{
  return idx_->get_gsite();
}

/*
const Field Dw::EvOd::md_force(const Field& eta,const Field& zeta) const{
  using namespace SUNmat_utils;

  int Nc = CommonPrms::instance()->Nc();
  int Nd = CommonPrms::instance()->Nd();

  Field et5 = gamma5(eta);
  Field zt5 = gamma5(zeta);

  Field fce(gf_->size());

  for(int mu=0; mu<Ndim_; ++mu){
    Field xie(fsize_), xz5(fsize_);
  
    sf_up_[mu]->setf(const_cast<Field&>(eta));
    (this->*mult_p[mu])(xie, sf_up_[mu]);

    sf_up_[mu]->setf(const_cast<Field&>(zt5));
    (this->*mult_p[mu])(xz5, sf_up_[mu]);

    for(int hs=0; hs<Nvol_; ++hs){
      SUNmat fe,fo;
      for(int a=0; a<Nc; ++a){
        for(int b=0; b<Nc; ++b){

          double fe_r = 0.0;
          double fe_i = 0.0;
          double fo_r = 0.0;
          double fo_i = 0.0;
	  
          for(int s=0; s<Nd; ++s){
	    size_t ra =ff_->index_r(a,s,hs);
	    size_t ia =ff_->index_i(a,s,hs);

	    size_t rb =ff_->index_r(b,s,hs);
	    size_t ib =ff_->index_i(b,s,hs);

	    fe_r += zeta[rb]*xie[ra] +zeta[ib]*xie[ia];
	    fe_i += zeta[rb]*xie[ia] -zeta[ib]*xie[ra];

	    fo_r += -xz5[rb]*et5[ra] -xz5[ib]*et5[ia];
	    fo_i += -xz5[rb]*et5[ia] +xz5[ib]*et5[ra];
          }
          fe.set(a,b,fe_r,fe_i);
          fo.set(a,b,fo_r,fo_i);
        }
      }
      fce.set(gf_->cslice(0,idx_->esec(hs),mu),anti_hermite(fe));
      fce.set(gf_->cslice(0,idx_->osec(hs),mu),anti_hermite(fo));
    }
  }
  fce *= -kpp_;
  return fce;
}  
*/

const vector<int> Dw::OdEv::get_gsite() const{
  return idx_->get_gsite();
}

const Field Dw::OdEv::mult(const Field& f)const{
  Field w(f.size());
  Dirac_Wilson::mult_core(w,const_cast<Field&>(f));
  return w;
}

const Field Dw::OdEv::mult_dag(const Field& f)const{
  return gamma5(mult(gamma5(f)));
}

const Field Dw::OdEv::md_force(const Field& eta,const Field& zeta) const{
  return Dirac_Wilson::md_force(eta,zeta);
}

/*
const Field Dw::OdEv::md_force(const Field& eta,const Field& zeta) const{
  using namespace SUNmat_utils;

  int Nc = CommonPrms::instance()->Nc();
  int Nd = CommonPrms::instance()->Nd();

  Field et5 = gamma5(eta);
  Field zt5 = gamma5(zeta);

  Field fce(gf_->size());

  for(int mu=0; mu<Ndim_; ++mu){
    Field xie(fsize_), xz5(fsize_);
  
    sf_up_[mu]->setf(const_cast<Field&>(eta));
    (this->*mult_p[mu])(xie, sf_up_[mu]);

    sf_up_[mu]->setf(const_cast<Field&>(zt5));
    (this->*mult_p[mu])(xz5, sf_up_[mu]);

    for(int hs=0; hs<Nvol_; ++hs){
      SUNmat fe,fo;
      for(int a=0; a<Nc; ++a){
        for(int b=0; b<Nc; ++b){

          double fe_r = 0.0;
          double fe_i = 0.0;
          double fo_r = 0.0;
          double fo_i = 0.0;
	  
          for(int s=0; s<Nd; ++s){
	    size_t ra =ff_->index_r(a,s,hs);
	    size_t ia =ff_->index_i(a,s,hs);

	    size_t rb =ff_->index_r(b,s,hs);
	    size_t ib =ff_->index_i(b,s,hs);

	    fe_r += -xz5[rb]*et5[ra] -xz5[ib]*et5[ia];
	    fe_i += -xz5[rb]*et5[ia] +xz5[ib]*et5[ra];

	    fo_r += zeta[rb]*xie[ra] +zeta[ib]*xie[ia];
	    fo_i += zeta[rb]*xie[ia] -zeta[ib]*xie[ra];
          }
          fe.set(a,b,fe_r,fe_i);
          fo.set(a,b,fo_r,fo_i);
        }
      }
      fce.set(gf_->cslice(0,idx_->esec(hs),mu),anti_hermite(fe));
      fce.set(gf_->cslice(0,idx_->osec(hs),mu),anti_hermite(fo));
    }
  }
  fce *= -kpp_;
  return fce;
}  
*/

const Field Dirac_Wilson_EvenOdd::mult(const Field& f) const{
  Field w(f);
  w -= Deo_.mult(Doe_.mult(f));
  return w;
}

const Field Dirac_Wilson_EvenOdd::mult_dag(const Field& f) const{
  Field w(f);
  w -= Deo_.mult_dag(Doe_.mult_dag(f));
  return w;
}

const Field Dirac_Wilson_EvenOdd::
md_force(const Field& eta,const Field& zeta) const{
  Field fce = Doe_.md_force(eta,Deo_.mult_dag(zeta));
  Field fce2 = Deo_.md_force(Doe_.mult(eta),zeta);
  CCIO::cout<<"Doe.fce: average_abs="<<fce.average_abs()<<std::endl;
  CCIO::cout<<"Deo.fce: average_abs="<<fce2.average_abs()<<std::endl;
  fce += fce2;
  return -fce;
}

const vector<int> Dirac_Wilson_EvenOdd::get_gsite() const {
  return Deo_.get_gsite();
}

