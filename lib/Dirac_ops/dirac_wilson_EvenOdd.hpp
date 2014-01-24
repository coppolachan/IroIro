/*! @file dirac_wilson_EvenOdd.hpp
  @brief Definition of Even Odd wilson operator
  Time-stamp: <2013-11-29 21:34:36 noaki>
*/
//---------------------------------------------------------
#ifndef DIRAC_WILSON_EVENODD_INCLUDED
#define DIRAC_WILSON_EVENODD_INCLUDED

#include "dirac_wilson.hpp"
#include "dirac_wilson_adjoint.hpp"

template<typename DWILSON>
class Dirac_Wilson_EvenOdd:public DiracWilsonLike_EvenOdd {
private:
  const DWILSON Deo_;
  const DWILSON Doe_;

  void md_force_eo(Field&,const Field&,const Field&)const;
  void md_force_oe(Field&,const Field&,const Field&)const;

  Dirac_Wilson_EvenOdd(const Dirac_Wilson_EvenOdd&);
  /*!< @brief simple copy is prohibited */
public:
  Dirac_Wilson_EvenOdd(double mass,const Field* u)
    :Deo_(mass,u,Dop::EOtag()),Doe_(mass,u,Dop::OEtag()){}

  Dirac_Wilson_EvenOdd(const XML::node& node,const Field* u)
    :Deo_(Dop::read_mass(node),u,Dop::EOtag()),
     Doe_(Dop::read_mass(node),u,Dop::OEtag()){}
  
  size_t fsize() const{return Deo_.fsize();}
  size_t gsize() const{return Deo_.gsize();}
  double getMass()const{return Deo_.getMass();}

  const Field gamma5(const Field& f) const{return Deo_.gamma5(f);}

  const Field mult(const Field&) const;
  const Field mult_dag(const Field&) const;

  const Field md_force(const Field&,const Field&) const;

  const Field mult_eo(const Field& f) const; 
  const Field mult_oe(const Field& f) const; 
  const Field mult_eo_dag(const Field& f) const;
  const Field mult_oe_dag(const Field& f) const;
  const Field mult_oo(const Field& f)const {return f;}
  const Field mult_ee(const Field& f)const {return f;}
  const Field mult_oo_inv(const Field& f)const {return f;}
  const Field mult_ee_inv(const Field& f)const {return f;}

  const DiracWilsonLike* getDeo()const{return &Deo_;}
  const DiracWilsonLike* getDoe()const{return &Doe_;}

  const Field* getGaugeField_ptr()const{ return Deo_.getGaugeField_ptr(); }
  
  double getKappa()const { return Deo_.getKappa();}
};

template <typename DWILSON>
const Field Dirac_Wilson_EvenOdd<DWILSON>::mult_eo(const Field& f) const{ 
  return Deo_.mult(f);
}

template <typename DWILSON>
const Field Dirac_Wilson_EvenOdd<DWILSON>::mult_oe(const Field& f) const{
  return Doe_.mult(f);
}

template <typename DWILSON>
const Field Dirac_Wilson_EvenOdd<DWILSON>::mult_eo_dag(const Field& f) const{
  return Doe_.mult_dag(f);
}

template <typename DWILSON>
const Field Dirac_Wilson_EvenOdd<DWILSON>::mult_oe_dag(const Field& f) const{
  return Deo_.mult_dag(f);
}

template <typename DWILSON>
const Field Dirac_Wilson_EvenOdd<DWILSON>::mult(const Field& f) const{
  Field w(f);
  w -= Deo_.mult(Doe_.mult(f));
  return w;
}

template <typename DWILSON>
const Field Dirac_Wilson_EvenOdd<DWILSON>::mult_dag(const Field& f) const{
  //return gamma5(mult(gamma5(f)));
  Field w(f);
  w -= mult_oe_dag(mult_eo_dag(f));
  return w;
}

template <typename DWILSON>
void Dirac_Wilson_EvenOdd<DWILSON>::
md_force_eo(Field& fce, const Field& eta,const Field& zeta) const{
  Deo_.md_force_p(fce,eta,zeta);
  Doe_.md_force_m(fce,eta,zeta);
}

template <typename DWILSON>
void Dirac_Wilson_EvenOdd<DWILSON>::
md_force_oe(Field& fce, const Field& eta,const Field& zeta) const{
  Doe_.md_force_p(fce,eta,zeta);
  Deo_.md_force_m(fce,eta,zeta);
}

template <typename DWILSON>
const Field Dirac_Wilson_EvenOdd<DWILSON>::
md_force(const Field& eta,const Field& zeta) const{
  Field fce(gsize());
  md_force_eo(fce,mult_oe(eta),zeta);
  md_force_oe(fce,eta,mult_eo_dag(zeta));
  fce *= getKappa();

  return fce;
}

#endif
