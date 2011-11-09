#include "dirac_wilson_eo.h"

const Dw::Even* Dirac_Wilson::get_even(){
  even_= new Dw::Even(kpp_,u_);
  return even_;
}

const Dw::Odd* Dirac_Wilson::get_odd(){
  odd_= new Dw::Odd(kpp_,u_);
  return odd_;
}

const Field Dw::Even::mult_oe(const Field& f) const{
  Field w(fsize_);
  mult_core(w,f);
  return w;
}
const Field Dw::Odd::mult_eo(const Field& f) const{
  Field w(fsize_);
  mult_core(w,f);
  return w;
}

