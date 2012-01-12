/*!
 * @file dirac_DomainWall_EvenOdd.hpp
 *
 * @brief Declaration of class Dirac_optimalDomainWall_EvenOdd (5d operator)
 */
#ifndef DIRAC_OPTIMALDOMAINWALL_EVENODD_INCLUDED
#define DIRAC_OPTIMALDOMAINWALL_EVENODD_INCLUDED

#include "dirac_wilson_DomainWall.h"
#include "dirac_wilson_EvenOdd.h"

/*!
 * @brief Defines the 5d Domain Wall operator with even/odd site indexing
 */
class Dirac_optimalDomainWall_EvenOdd : public DiracWilsonLike_EvenOdd {
  const Dirac_Wilson_EvenOdd* Dw_;/*!< @brief Dirac Kernel - Wilson operator */ 
  Dirac_optimalDomainWall_params Params;
  Preconditioner* Precond_;

  size_t N5_;/*!< @brief Length of 5th dimension */
  size_t f4size_;
  size_t fsize_;
  size_t gsize_;
  const double M0_;

  //declaration of concrete preconditioners
  class NoPrecond: public Preconditioner {
    Dirac_optimalDomainWall_EvenOdd* DWF_;
  public: 
    NoPrecond(Dirac_optimalDomainWall_EvenOdd* DWF): DWF_(DWF){}
    const Field precondition(const Field&) const;  
    const Field mult(const Field&) const;  
    const Field mult_dag(const Field&) const;
    const Field left(const Field&) const;  
    const Field right(const Field&) const;  
  };
  
  class LUPrecond : public Preconditioner {
    Dirac_optimalDomainWall_EvenOdd* DWF_;
  public: 
    LUPrecond(Dirac_optimalDomainWall_EvenOdd* DWF): DWF_(DWF){}
    const Field precondition(const Field&) const;  
    const Field mult(const Field&) const;  
    const Field mult_dag(const Field&) const;
    const Field left(const Field&) const;  
    const Field right(const Field&) const; 
  };

  Preconditioner* choose_Preconditioner(int PrecondID);

  const Field get4d(const Field& f5,int s) const{
    return Field(f5[std::slice(s*f4size_,f4size_,1)]);
  }
  void set5d(Field& f5,const Field& f4,int s) const{
    f5.set(std::slice(s*f4size_,f4size_,1),f4.getva());
  }
  void add5d(Field& f5,const Field& f4,int s) const{
    f5.add(std::slice(s*f4size_,f4size_,1),f4.getva());
  }

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
    :Params(Dcopy.Params, Type),
     Dw_(Dcopy.Dw_),
     N5_(Dcopy.N5_),
     f4size_(Dcopy.f4size_),
     fsize_(Dcopy.fsize_),
     gsize_(Dcopy.gsize_),
     M0_(Dcopy.M0_),
     Precond_(choose_Preconditioner(Params.Preconditioning_))//cannot just copy
  {}

  ~Dirac_optimalDomainWall_EvenOdd(){
    #if VERBOSITY>4
    CCIO::cout << "Deleting Dirac_optimalDomainWall" << std::endl;
    #endif
    
    delete Precond_; 
  }
  
  size_t f4size() const{ return f4size_;}
  size_t fsize()  const{ 
    CCIO::cout<<"Dirac_optimalDomainWall::fsize_="
	      <<fsize_<<std::endl;
    return fsize_; }
  size_t gsize()  const{ return gsize_; }
  
  const Field operator()(int, const Field&) const{}

  const double getMass() const{return Params.mq_;}

  const Field gamma5_4d(const Field& f4d) const{return Dw_->gamma5(f4d);}
  
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
