/*! @file dirac_wilson_adjoint_EvenOdd.hpp
  @brief Definition of Even Odd wilson operator with the adjoint rep.
  Time-stamp: <2013-11-29 16:11:08 noaki>
*/
//---------------------------------------------------------
#ifndef DIRAC_WILSON_ADJOINT_EVENODD_INCLUDED
#define DIRAC_WILSON_ADJOINT_EVENODD_INCLUDED

#include "dirac_wilson_adjoint.hpp"

class Dirac_Wilson_Adjoint_EvenOdd :public DiracWilsonLike_EvenOdd {
private:
  const Dirac_Wilson_Adjoint Deo_;
  const Dirac_Wilson_Adjoint Doe_;

  void md_force_eo(Field&,const Field&,const Field&)const;
  void md_force_oe(Field&,const Field&,const Field&)const;

  Dirac_Wilson_Adjoint_EvenOdd(const Dirac_Wilson_Adjoint_EvenOdd&);
  /*!< @brief simple copy is prohibited */
public:
  Dirac_Wilson_Adjoint_EvenOdd(double mass,const Field* u)
    :Deo_(mass,u,Dop::EOtag()),Doe_(mass,u,Dop::OEtag()){}

  Dirac_Wilson_Adjoint_EvenOdd(const XML::node& node,const Field* u)
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

#endif
