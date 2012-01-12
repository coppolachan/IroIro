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

  Dirac_optimalDomainWall_params Params;
  Preconditioner* Precond_;
  const double M0_;

  Preconditioner* choose_Preconditioner(int PrecondID);

public:
  Dirac_optimalDomainWall_EvenOdd(XML::node DWF_node,
				  const Dirac_Wilson_EvenOdd* Kernel);
  
  Dirac_optimalDomainWall_EvenOdd(Dirac_optimalDomainWall_params Prm,
				  const Dirac_Wilson_EvenOdd* Kernel);
 
  Dirac_optimalDomainWall_EvenOdd(const double b,
				  const double c,
				  const double mq,
				  const std::vector<double>& omega,
				  const Dirac_Wilson_EvenOdd* Kernel,
				  Preconditioners Precond = NoPreconditioner);

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
  const Field gamma5_4d(const Field& f4d) const{return Deo_->gamma5(f4d);}
  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;

  //Preconditioning methods
  const Field mult_prec    (const Field& in)const{return Precond_->mult(in);}
  const Field mult_dag_prec(const Field& in)const{return Precond_->mult_dag(in);}
  const Field left_precond (const Field& in)const{return Precond_->left(in);}
  const Field right_precond(const Field& in)const{return Precond_->right(in);}
  //////////////////////////////////////////////////////////////////////

  const Field Dminus(const Field&)const;
  const Field gamma5(const Field&) const;

  /*!
   * @brief Calculates the \f$L_+(m)\f$
   */
  const Field proj_p(const Field& f4) const;
    
  /*!
   * @brief Calculates the \f$L_-(m)\f$
   */
  const Field proj_m(const Field& f4) const;

  const Field Bproj(const Field& v5d) const;
  const Field Bproj_dag(const Field& v4d) const;
  const Field R5(const Field&) const;
  const Field R5g5(const Field&) const;

  const Field md_force( const Field& eta,const Field& zeta) const;
  const Format::Format_F get_fermionFormat() const{
    Format::Format_F ff = Dw_->get_fermionFormat();
    return Format::Format_F(ff.Nvol(),N5_);
  }
  const std::vector<int> get_gsite() const { return Dw_->get_gsite();}
};


#endif
