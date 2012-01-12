//--------------------------------------------------------------
/*
  @file dirac_wilson_EvenOdd.cpp
  @brief Definition of Even Odd wilson operator
*/
//--------------------------------------------------------------
#include "dirac_wilson_EvenOdd.h"
using namespace std;

const Field Dirac_Wilson_EvenOdd::mult_eo(const Field& f) const{ 
  Field w(f.size());
  Deo_.mult_core(w,const_cast<Field&>(f));
  return w;
}

const Field Dirac_Wilson_EvenOdd::mult_oe(const Field& f) const{
  Field w(f.size());
  Doe_.mult_core(w,const_cast<Field&>(f));
  return w;
}

const Field Dirac_Wilson_EvenOdd::mult_eo_dag(const Field& f) const{
  return gamma5(mult_oe(gamma5(f)));
}

const Field Dirac_Wilson_EvenOdd::mult_oe_dag(const Field& f) const{
  return gamma5(mult_eo(gamma5(f)));
}

const Field Dirac_Wilson_EvenOdd::mult(const Field& f) const{
  Field w(f);
  w -= mult_eo(mult_oe(f));
  return w;
}

const Field Dirac_Wilson_EvenOdd::mult_dag(const Field& f) const{
  //return gamma5(mult(gamma5(f)));
  Field w(f);
  w -= mult_oe_dag(mult_eo_dag(f));
  return w;
}

const vector<int> Dirac_Wilson_EvenOdd::get_gsite() const {
  return SiteIndex_eo::instance()->get_gsite();
}

const Field Dirac_Wilson_EvenOdd::
md_force(const Field& eta,const Field& zeta) const{

  Field fce(gsize());
  Deo_.md_force_p(fce,mult_oe(eta),zeta);
  Doe_.md_force_m(fce,mult_oe(eta),zeta);
  Deo_.md_force_m(fce,eta,mult_eo_dag(zeta));
  Doe_.md_force_p(fce,eta,mult_eo_dag(zeta));
  
  fce *= getKappa();
  return fce;
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
  fe.set(osec,mult_oe(eta).getva());
  fz.set(osec,mult_eo_dag(zeta).getva());

  Field fce(Dw_.md_force(fe,fz));
  return -fce;
}
*/


