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

void Dw::EvOd::md_force_p(Field& fce,const Field& eta,const Field& zeta)const{
  Dirac_Wilson::md_force_p(fce,eta,zeta);
}
void Dw::EvOd::md_force_m(Field& fce,const Field& eta,const Field& zeta)const{
  Dirac_Wilson::md_force_m(fce,eta,zeta);
}

const vector<int> Dw::EvOd::get_gsite() const{ return idx_->get_gsite();}

const Field Dw::OdEv::mult(const Field& f)const{
  Field w(f.size());
  Dirac_Wilson::mult_core(w,const_cast<Field&>(f));
  return w;
}

const Field Dw::OdEv::mult_dag(const Field& f)const{
  return gamma5(mult(gamma5(f)));
}

void Dw::OdEv::md_force_p(Field& fce,const Field& eta,const Field& zeta)const{
  Dirac_Wilson::md_force_p(fce,eta,zeta);
}
void Dw::OdEv::md_force_m(Field& fce,const Field& eta,const Field& zeta)const{
  Dirac_Wilson::md_force_m(fce,eta,zeta);
}

const vector<int> Dw::OdEv::get_gsite() const{return idx_->get_gsite();}

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

/*
const Field Dirac_Wilson_EvenOdd::
md_force(const Field& eta,const Field& zeta) const{

  Format::Format_F ff(fsize()*2);
  Field fe(ff.size());
  Field fz(ff.size());

  SiteIndex_eo* idx = SiteIndex_eo::instance();
  
  valarray<size_t> esec(ff.get_sub(idx->esec()));
  fe.set(esec,eta.getva());
  fz.set(esec,zeta.getva());

  valarray<size_t> osec(ff.get_sub(idx->osec()));
  fe.set(osec,Doe_.mult(eta).getva());
  fz.set(osec,Deo_.mult_dag(zeta).getva());

  Field fce(Dw_.md_force(fe,fz));
  return -fce;
}
*/

const Field Dirac_Wilson_EvenOdd::
md_force(const Field& eta,const Field& zeta) const{

  Field fce(gsize());
  Deo_.md_force_p(fce,Doe_.mult(eta),zeta);
  Doe_.md_force_m(fce,Doe_.mult(eta),zeta);
  Deo_.md_force_m(fce,eta,Deo_.mult_dag(zeta));
  Doe_.md_force_p(fce,eta,Deo_.mult_dag(zeta));
  
  fce *= getKappa();
  return fce;
}

