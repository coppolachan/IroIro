//--------------------------------------------------------------
/*
  @file dirac_wilson_EvenOdd.cpp
  @brief Definition of Even Odd wilson operator
*/
//--------------------------------------------------------------
#include "dirac_wilson_EvenOdd.h"

const Field Dw::EvOd::mult(const Field& f)const{
  Field w(f.size());
  Dirac_Wilson::mult_core(w,const_cast<Field&>(f));
  return w;
}

const Field Dw::EvOd::mult_dag(const Field& f)const{
  return gamma5(mult(gamma5(f)));
}

const Field Dw::OdEv::mult(const Field& f)const{
  Field w(f.size());
  Dirac_Wilson::mult_core(w,const_cast<Field&>(f));
  return w;
}

const Field Dw::OdEv::mult_dag(const Field& f)const{
  return gamma5(mult(gamma5(f)));
}

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
  Field fce = Deo_.md_force(eta,Doe_.mult(zeta));
  fce += Doe_.md_force(Deo_.mult_dag(eta),zeta);
  return -fce;
}

