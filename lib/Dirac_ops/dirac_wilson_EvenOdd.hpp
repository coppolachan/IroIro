//---------------------------------------------------------
/*
  @file dirac_wilson_EvenOdd.hpp
  @brief Definition of Even Odd wilson operator
*/
//---------------------------------------------------------
#ifndef DIRAC_WILSON_EVENODD_INCLUDED
#define DIRAC_WILSON_EVENODD_INCLUDED

#include "dirac_wilson.hpp"

class Dirac_Wilson_EvenOdd:public DiracWilsonLike_EvenOdd {
private:
  const Dirac_Wilson Deo_;
  const Dirac_Wilson Doe_;
  //const Dirac_Wilson Dw_;

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

  const Field gamma5(const Field& f) const{return Deo_.gamma5(f);}

  const Field mult(const Field&) const;
  const Field mult_dag(const Field&) const;

  ////////////////////////////////////////Preconditioned versions
  // EvenOdd operator has no preconditioner now 
  const Field mult_prec     (const Field&f)const{return f;}
  const Field mult_dag_prec (const Field&f)const{return f;}
  const Field left_prec     (const Field&f)const{return f;}
  const Field right_prec    (const Field&f)const{return f;}
  const Field left_dag_prec (const Field&f)const{return f;}
  const Field right_dag_prec(const Field&f)const{return f;}
  //////////////////////////////////////////////////////////////
  
  const Field md_force(const Field&,const Field&) const;
  void update_internal_state(){}

  const Field mult_eo(const Field& f) const; 
  const Field mult_oe(const Field& f) const; 
  const Field mult_eo_dag(const Field& f) const;
  const Field mult_oe_dag(const Field& f) const;
  const Field mult_oo(const Field& f)const {return f;}
  const Field mult_ee(const Field& f)const {return f;}
  const Field mult_oo_inv(const Field& f)const {return f;}
  const Field mult_ee_inv(const Field& f)const {return f;}

  void get_RandGauss(std::valarray<double>& phi,const RandNum& rng)const;
  double getKappa()const { return Deo_.getKappa();}
  int Nvol()const{ return Deo_.Nvol();}
};

#endif
