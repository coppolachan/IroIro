/*!
 * @file dirac_DomainWall_EvenOdd.hpp
 *
 * @brief Declaration of class Dirac_optimalDomainWall_EvenOdd (5d operator)
 */
#ifndef DIRAC_OPTIMALDOMAINWALL_EVENODD_INCLUDED
#define DIRAC_OPTIMALDOMAINWALL_EVENODD_INCLUDED

#include "dirac_wilson_DomainWall.h"

/*!
 * @brief Defines the 5d Domain Wall operator with even/odd site indexing
 */
class Dirac_optimalDomainWall_EvenOdd : public DiracWilsonLike_EvenOdd {
  const Dirac_optimalDomainWall Deo_;
  const Dirac_optimalDomainWall Doe_;

public:
  Dirac_optimalDomainWall_EvenOdd(XML::node DWF_node,
				  const Dirac_Wilson* Kernel_eo,
				  const Dirac_Wilson* Kernel_oe)
  :Deo_(DWF_node,Kernel_eo),Doe_(DWF_node,Kernel_oe){}
  
  Dirac_optimalDomainWall_EvenOdd(const double b,
				  const double c,
				  const double mq,
				  const std::vector<double>& omega,
				  const Dirac_Wilson* Kernel_eo,
				  const Dirac_Wilson* Kernel_oe,
				  Preconditioners Precond = NoPreconditioner)
    :Deo_(b,c,mq,omega,Kernel_eo,Precond),
     Doe_(b,c,mq,omega,Kernel_oe,Precond){}

  /*! @brief Copy constructor to build the Pauli-Villars operator */
  Dirac_optimalDomainWall_EvenOdd(const Dirac_optimalDomainWall_EvenOdd& Dcopy, 
				  DWFType Type = Standard)
    :Deo_(Dcopy.Deo_,Type),Doe_(Dcopy.Doe_,Type){}
  
  ~Dirac_optimalDomainWall_EvenOdd(){}
  
  size_t f4size()const{ return Deo_.f4size();}
  size_t fsize() const{ return Deo_.fsize(); }
  size_t gsize() const{ return Deo_.gsize(); }
  
  const Field operator()(int, const Field&) const{}
  const double getMass() const{return Deo_.getMass();}
  const Field gamma5_4d(const Field& f4d) const{return Deo_.gamma5(f4d);}
  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;

  //Preconditioning methods
  const Field mult_prec    (const Field& f)const{return f;}
  const Field mult_dag_prec(const Field& f)const{return f;}
  const Field left_precond (const Field& f)const{return f;}
  const Field right_precond(const Field& f)const{return f;}
  //////////////////////////////////////////////////////////////////////

  const Field gamma5(const Field&) const;
  const Field md_force( const Field& eta,const Field& zeta) const;

  const Field mult_eo(const Field& f) const; 
  const Field mult_oe(const Field& f) const; 
  const Field mult_eo_dag(const Field& f) const;
  const Field mult_oe_dag(const Field& f) const;
  const Field mult_oo(const Field& f)const;
  const Field mult_ee(const Field& f)const;
  const Field mult_oo_inv(const Field& f)const;
  const Field mult_ee_inv(const Field& f)const;

  const Format::Format_F get_fermionFormat() const{
    Format::Format_F ff = Deo_->get_fermionFormat();
    return Format::Format_F(ff.Nvol(),N5_);
  }
  const std::vector<int> get_gsite() const { return Deo_.get_gsite();}
};


#endif
