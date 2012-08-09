//------------------------------------------------------------------
/*!@file   dirac_staggered_EvenOdd.hpp
 * @brief  Definition of the staggered operator on even/odd sites
 */
//------------------------------------------------------------------
#ifndef DIRAC_STAGGERED_EVENODD_INCLUDED
#define DIRAC_STAGGERED_EVENODD_INCLUDED

#include "dirac_wilson_EvenOdd.hpp"

class Dirac_staggered_EvenOdd :public DiracStaggeredLike_EvenOdd{
private:
  const Dirac_staggered Deo_;
  const Dirac_staggered Doe_;
  
  void md_force_eo(Field&,const Field&,const Field&)const;
  void md_force_oe(Field&,const Field&,const Field&)const;

  Dirac_staggered_EvenOdd(const Dirac_staggered_EvenOdd&);
  /*!< @brief simple copy is prohibited */
public:
  // manual construction
  Dirac_staggered_EvenOdd(double mass,const Field* u)
    :Deo_(mass,u,Dop::EOtag()),Doe_(mass,u,Dop:OEtag()){}

  // xml construction 
  Dirac_Wilson_EvenOdd(const XML::node& node,const Field* u)
    :Deo_(Dop::read_mass(node),u,Dop::EOtag()),
     Doe_(Dop::read_mass(node),u,Dop::OEtag()){}
  
  size_t fsize() const{return Deo_.fsize();}
  size_t gsize() const{return Deo_.gsize();}

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

  double get_mq()const{ return Deo_.get_mq();}
  const std::vector<int> get_gsite() const;
};
#endif
