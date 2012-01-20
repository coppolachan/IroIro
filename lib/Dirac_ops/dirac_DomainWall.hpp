/*!
 * @file dirac_DomainWall.hpp
 *
 * @brief Declaration of class Dirac_optimalDomainWall (5d operator)
 */
#ifndef DIRAC_OPTIMALDOMAINWALL_INCLUDED
#define DIRAC_OPTIMALDOMAINWALL_INCLUDED

#include <string>
#include <string.h>
#include "include/field.h"
#include "include/macros.hpp"

#include "include/pugi_interface.h"
#include "Tools/EnumToString.hpp"

#include "dirac_wilson.h"
#include "dirac_Preconditioners.hpp"

namespace DomainWallFermions {

  struct EvenOdd_tag{};

  const std::vector<double> 
  getOmega(int Ns,double lambda_min,double lambda_max);
}

enum DWFType {Standard, PauliVillars};
enum Preconditioners {NoPreconditioner, LUPreconditioner};
Begin_Enum_String( Preconditioners ) {
  Enum_String( NoPreconditioner );
  Enum_String( LUPreconditioner );
} 
End_Enum_String;

/*!
 * @brief Container for parameter of the 5d Optimal Domain Wall operator
 * Parameters for \f[D_{\rm dwf}(m_q) =\omega D_W(-m_0)(1+cL(m_q))+(1-L(m_q))\f]
 */
struct Dirac_optimalDomainWall_params{
  double N5dim_;/*!< @brief the length in the 5th direction (must be even) */
  double b_;    /*!< @brief scale factor (b!=1 for scaled Shamir H_T) */
  double c_;    /*!< @brief the kernel (H_W (c=0) or H_T (c=1)) */
  double mq_;   /*!< @brief Bare quark mass \f$m_q\f$ */
  Preconditioners Preconditioning_; /*!< @brief Name of the preconditioner */
  std::vector<double> omega_;/*!< @brief Weights defining the approximation */
  std::vector<double> bs_, cs_;
  std::vector<double> dp_, dm_;
  std::vector<double> es_, fs_;

  Dirac_optimalDomainWall_params(XML::node DWF_node);

  Dirac_optimalDomainWall_params(const double b,
				 const double c,
				 const double mass,
				 const std::vector<double> omega);

  Dirac_optimalDomainWall_params(const Dirac_optimalDomainWall_params& Par,
				 DWFType Type = Standard);
};

/*!
 * @brief Defines the 5d Domain Wall operator
 *
 * Defines \f[D_{\rm dwf}(m_q) = \omega D_W(-m_0)(1+cL(m_q))+(1-L(m_q))\f]
 * with:
 * \f[L(m) = P_+ L_+(m) + P_- L_-(m)\f]
 * using proj_p() e proj_m() methods
 *
 * Kernel is given by
 * \f[H = \frac{\gamma_5 b D_W}{ 2 + c D_W } \f]
 *
 */
class Dirac_optimalDomainWall : public DiracWilsonLike {
  const Dirac_Wilson* Dw_; /*!< @brief Dirac Kernel - Wilson operator */ 
  Dirac_optimalDomainWall_params Params;
  Preconditioner* Precond_;

  size_t N5_;/*!< @brief Length of 5th dimension */
  size_t f4size_;
  size_t fsize_;
  size_t gsize_;
  const double M0_;

  //declaration of concrete preconditioners
  class NoPrecond: public Preconditioner {
    Dirac_optimalDomainWall* DWF_;
  public: 
    NoPrecond(Dirac_optimalDomainWall* DWF): DWF_(DWF){}
    const Field mult(const Field& f5) const{return DWF_->mult(f5);}
    const Field mult_dag(const Field& f5) const{return DWF_->mult_dag(f5);}
    const Field left(const Field& f5) const{ return f5;}
    const Field right(const Field& f5) const{ return f5;}
    const Field left_dag(const Field&) const;  
    const Field right_dag(const Field&) const;  
  };
  
  class LUPrecond : public Preconditioner {
    Dirac_optimalDomainWall* DWF_;
    const Field LU(const Field& f5) const{return DWF_->mult_hq(f5);}
    const Field LU_inv(const Field& f5) const{return DWF_->mult_hq_inv(f5);}
    const Field LU_dag(const Field& f5) const{return DWF_->mult_hq_dag(f5);}
    const Field LU_dinv(const Field& f5)const{return DWF_->mult_hq_dinv(f5);}
  public: 
    LUPrecond(Dirac_optimalDomainWall* DWF): DWF_(DWF){}
    const Field mult(const Field&) const;  
    const Field mult_dag(const Field&) const;
    const Field left(const Field& f5) const{return LU_inv(f5);}
    const Field right(const Field& f5) const{ return f5;}
    const Field left_dag(const Field&) const;  
    const Field right_dag(const Field&) const; 
  };

  void set_arrays();
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

  void mult_a0(Field&,const Field&)const;
  void mult_a1(Field&,const Field&)const;

  void mult_dag_a0(Field&,const Field&)const;
  void mult_dag_a1(Field&,const Field&)const;

  void(Dirac_optimalDomainWall::*mult_core)(Field&,const Field&)const;
  void(Dirac_optimalDomainWall::*mult_dag_core)(Field&,const Field&)const;

public:

  /*! @brief constructors to create an instance with normal indexing */
  Dirac_optimalDomainWall(XML::node DWF_node,
			  const Dirac_Wilson* Kernel)
    :Params(Dirac_optimalDomainWall_params(DWF_node)),
     Dw_(Kernel),
     N5_(Params.omega_.size()),
     f4size_(Dw_->fsize()),
     fsize_(f4size_*N5_),
     gsize_(Dw_->gsize()),
     M0_(1.0/(2.0*(Dw_->getKappa()))-4.0),
     mult_core(&Dirac_optimalDomainWall::mult_a1),
     mult_dag_core(&Dirac_optimalDomainWall::mult_dag_a1),
     Precond_(choose_Preconditioner(Params.Preconditioning_)){
    set_arrays();
    CCIO::cout << "Dirac_optimalDomainWall::N5_="<<N5_<<std::endl;
#if VERBOSITY>4
    CCIO::cout << "Created Dirac_optimalDomainWall" << std::endl;
#endif
  }
  
  Dirac_optimalDomainWall(Dirac_optimalDomainWall_params Prm,
			  const Dirac_Wilson* Kernel)
    :Params(Prm),
     Dw_(Kernel),
     N5_(Params.omega_.size()),
     f4size_(Dw_->fsize()),
     fsize_(f4size_*N5_),
     gsize_(Dw_->gsize()),
     M0_(1.0/(2.0*(Dw_->getKappa()))-4.0),
     mult_core(&Dirac_optimalDomainWall::mult_a1),
     mult_dag_core(&Dirac_optimalDomainWall::mult_dag_a1),
     Precond_(choose_Preconditioner(Params.Preconditioning_)){
    set_arrays();
  } 
 
  Dirac_optimalDomainWall(const double b,
			  const double c,
			  const double mq,
			  const std::vector<double>& omega,
			  const Dirac_Wilson* Kernel,
			  Preconditioners Precond = NoPreconditioner)
    :Params(b,c,mq,omega),
     Dw_(Kernel),
     N5_(Params.omega_.size()),
     f4size_(Dw_->fsize()),
     fsize_(f4size_*N5_),
     gsize_(Dw_->gsize()),
     M0_(1.0/(2.0*(Dw_->getKappa()))-4.0),
     mult_core(&Dirac_optimalDomainWall::mult_a1),
     mult_dag_core(&Dirac_optimalDomainWall::mult_dag_a1),
     Precond_(choose_Preconditioner(Precond)){
    set_arrays();
  }

  /*! @brief constructors to create an instance with the e/o indexing */
  Dirac_optimalDomainWall(XML::node DWF_node,
			  const Dirac_Wilson* Kernel,
			  DomainWallFermions::EvenOdd_tag)
    :Params(Dirac_optimalDomainWall_params(DWF_node)),
     Dw_(Kernel),
     N5_(Params.omega_.size()),
     f4size_(Dw_->fsize()),
     fsize_(f4size_*N5_),
     gsize_(Dw_->gsize()),
     M0_(1.0/(2.0*(Dw_->getKappa()))-4.0),
     mult_core(&Dirac_optimalDomainWall::mult_a0),
     mult_dag_core(&Dirac_optimalDomainWall::mult_dag_a0),
     Precond_(NULL){
    set_arrays();
  }

  Dirac_optimalDomainWall(Dirac_optimalDomainWall_params Prm,
			  const Dirac_Wilson* Kernel,
			  DomainWallFermions::EvenOdd_tag)
    :Params(Prm),
     Dw_(Kernel),
     N5_(Params.omega_.size()),
     f4size_(Dw_->fsize()),
     fsize_(f4size_*N5_),
     gsize_(Dw_->gsize()),
     M0_(1.0/(2.0*(Dw_->getKappa()))-4.0),
     mult_core(&Dirac_optimalDomainWall::mult_a0),
     mult_dag_core(&Dirac_optimalDomainWall::mult_dag_a0),
     Precond_(NULL){
    set_arrays();
  } 
 
  Dirac_optimalDomainWall(const double b,
			  const double c,
			  const double mq,
			  const std::vector<double>& omega,
			  const Dirac_Wilson* Kernel,
			  DomainWallFermions::EvenOdd_tag)
    :Params(b,c,mq,omega),
     Dw_(Kernel),
     N5_(Params.omega_.size()),
     f4size_(Dw_->fsize()),
     fsize_(f4size_*N5_),
     gsize_(Dw_->gsize()),
     M0_(1.0/(2.0*(Dw_->getKappa()))-4.0),
     mult_core(&Dirac_optimalDomainWall::mult_a0),
     mult_dag_core(&Dirac_optimalDomainWall::mult_dag_a0),
     Precond_(NULL){
    set_arrays();
  }

  /*! @brief Copy constructor to build the Pauli-Villars operator */
  Dirac_optimalDomainWall(const Dirac_optimalDomainWall& Dcopy, 
			  DWFType Type = Standard)
    :Params(Dcopy.Params, Type),
     Dw_(Dcopy.Dw_),
     N5_(Dcopy.N5_),
     f4size_(Dcopy.f4size_),
     fsize_(Dcopy.fsize_),
     gsize_(Dcopy.gsize_),
     M0_(Dcopy.M0_),
     mult_core(&Dirac_optimalDomainWall::mult_a1),
     mult_dag_core(&Dirac_optimalDomainWall::mult_dag_a1),
     Precond_(choose_Preconditioner(Params.Preconditioning_)){}

  Dirac_optimalDomainWall(const Dirac_optimalDomainWall& Dcopy, 
			  DomainWallFermions::EvenOdd_tag,
			  DWFType Type = Standard)
    :Params(Dcopy.Params, Type),
     Dw_(Dcopy.Dw_),
     N5_(Dcopy.N5_),
     f4size_(Dcopy.f4size_),
     fsize_(Dcopy.fsize_),
     gsize_(Dcopy.gsize_),
     M0_(Dcopy.M0_),
     mult_core(&Dirac_optimalDomainWall::mult_a0),
     mult_dag_core(&Dirac_optimalDomainWall::mult_dag_a0),
     Precond_(NULL){}
  
  //////////////////////////////////////////
  ~Dirac_optimalDomainWall(){
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

  // mult in the heavy quark limit
  const Field mult_hq(const Field& f5) const;
  const Field mult_hq_inv(const Field& f5) const;
  const Field mult_hq_dag(const Field& f5) const;
  const Field mult_hq_dinv(const Field& f5) const;

  //Preconditioning methods
  const Field mult_prec    (const Field& in)const{return Precond_->mult(in);}
  const Field mult_dag_prec(const Field& in)const{return Precond_->mult_dag(in);}
  const Field left_prec (const Field& in)const{return Precond_->left(in);}
  const Field right_prec(const Field& in)const{return Precond_->right(in);}
  const Field left_dag_prec (const Field& in)const{
    return Precond_->left_dag(in);}
  const Field right_dag_prec(const Field& in)const{
    return Precond_->right_dag(in);}
  //////////////////////////////////////////////////////////////////////

  const Field Dminus(const Field&)const;
  const Field gamma5(const Field&) const;

  /*! @brief Calculates the \f$L_+(m)\f$ */
  const Field proj_p(const Field& f4) const;
    
  /*! @brief Calculates the \f$L_-(m)\f$ */
  const Field proj_m(const Field& f4) const;

  const Field Bproj(const Field& v5d) const;
  const Field Bproj_dag(const Field& v4d) const;
  const Field R5(const Field&) const;
  const Field R5g5(const Field&) const;

  void md_force_p(Field&,const Field&,const Field&) const;
  void md_force_m(Field&,const Field&,const Field&) const;

  const Field md_force( const Field& eta,const Field& zeta) const;
  const Format::Format_F get_fermionFormat() const{
    Format::Format_F ff = Dw_->get_fermionFormat();
    return Format::Format_F(ff.Nvol(),N5_);
  }
  const std::vector<int> get_gsite() const { return Dw_->get_gsite();}


  void update_internal_state(){};

};


#endif
