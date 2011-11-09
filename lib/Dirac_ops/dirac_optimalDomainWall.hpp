/*!
 * @file dirac_optimalDomainWall.hpp
 *
 * @brief Declaration of class Dirac_optimalDomainWall (5d operator)
 */
#ifndef DIRAC_OPTIMALDOMAINWALL_INCLUDED
#define DIRAC_OPTIMALDOMAINWALL_INCLUDED

#include "include/field.h"
#include "include/pugi_interface.h"
#include "Dirac_ops/dirac_wilson.h"


enum {Standard, PauliVillars};

/*!
 * @brief Container for parameter of the 5d Optimal Domain Wall operator
 *
 * Parameters for \f[D_{\rm dwf}(m_q) = \omega D_W(-m_0)(1+cL(m_q))+(1-L(m_q))\f]
 * 
 *
 */
struct Dirac_optimalDomainWall_params{
  double c_;/*!< @brief With omega_ selects the type of domain wall fermion
	     *
	     * For \f$c = 0, \omega_s = 1\f$ conventional 
	     * domain-wall fermion is chosen
	     */
  double mq_;/*!< @brief Bare quark mass \f$m_q\f$ */
  std::vector<double> omega_;/*!< @brief Diagonal matrix
			      * \f$\omega = \{\omega_s\}\f$
			      */

  Dirac_optimalDomainWall_params(XML::node DWF_node){
    double omega_diag;
    XML::read(DWF_node, "c", c_);
    XML::read(DWF_node, "mass", mq_);
    XML::read(DWF_node, "omega", omega_diag);
    //fill the vector omega_ with omega_diag
  }

  Dirac_optimalDomainWall_params(const double c,
				 const double mass,
				 const std::vector<double> omega){
    c_ = c;
    mq_ = mass;
    omega_ = omega;
  }

  Dirac_optimalDomainWall_params(const Dirac_optimalDomainWall_params& Par,
				 int Type = 0){
    switch (Type){
    case Standard:
      c_ = Par.c_;
      mq_ = Par.mq_;
      omega_ = Par.omega_;
      break;
    case PauliVillars:
      c_ = Par.c_;
      mq_ = 1.0;
      omega_ = Par.omega_;
      break;
    default:
      abort();
    }
  }



};


/*!
 * @brief Defines the 5d Optimal Domain Wall operator
 *
 * Defines \f[D_{\rm dwf}(m_q) = \omega D_W(-m_0)(1+cL(m_q))+(1-L(m_q))\f]
 * with:
 * \f[L(m) = P_+ L_+(m) + P_- L_-(m)\f]
 * using proj_p() e proj_m() methods
 *
 */
class Dirac_optimalDomainWall : public DiracWilsonLike {
  const Dirac_Wilson* Dw_; /*!< Dirac Kernel - Wilson operator */ 
  const Dirac_optimalDomainWall_params Params;

  size_t N5_;/*!< Length of 5th dimension */
  size_t f4size_;
  size_t fsize_;
  size_t gsize_;
  const double M0_;

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
  Dirac_optimalDomainWall(XML::node DWF_node,
			  const Dirac_Wilson* Kernel)
    :Params(Dirac_optimalDomainWall_params(DWF_node)),
     Dw_(Kernel),
     N5_(Params.omega_.size()),
     f4size_(Dw_->fsize()),
     fsize_(f4size_*N5_),
     gsize_(Dw_->gsize()),
     M0_(1.0/(2.0*(Dw_->getKappa()))-4.0){
  }
  
  Dirac_optimalDomainWall(Dirac_optimalDomainWall_params Prm,
			  const Dirac_Wilson* Kernel)
    :Params(Prm),
     Dw_(Kernel),
     N5_(Params.omega_.size()),
     f4size_(Dw_->fsize()),
     fsize_(f4size_*N5_),
     gsize_(Dw_->gsize()),
     M0_(1.0/(2.0*(Dw_->getKappa()))-4.0){} 

 
  Dirac_optimalDomainWall(const double c,
			  const double mq,
			  const std::vector<double>& omega,
			  const Dirac_Wilson* Kernel)
    :Params(c,mq,omega),
     Dw_(Kernel),
     N5_(Params.omega_.size()),
     f4size_(Dw_->fsize()),
     fsize_(f4size_*N5_),
     gsize_(Dw_->gsize()),
     M0_(1.0/(2.0*(Dw_->getKappa()))-4.0){}


  /*! @brief Copy constructor to build the Pauli-Villars operator */
  Dirac_optimalDomainWall(const Dirac_optimalDomainWall& Dcopy, 
			  int Type = 0)
    :Params(Dcopy.Params, Type),
     Dw_(Dcopy.Dw_),
     N5_(Dcopy.N5_),
     f4size_(Dcopy.f4size_),
     fsize_(Dcopy.fsize_),
     gsize_(Dcopy.gsize_),
     M0_(Dcopy.M0_){}

  ~Dirac_optimalDomainWall(){}
  
  size_t f4size() const{ return f4size_;}
  size_t fsize()  const{ return fsize_; }
  size_t gsize()  const{ return gsize_; }
  
  const Field operator()(int, const Field&) const{};

  const double getMass() const{return Params.mq_;}

  const Field gamma5_4d(const Field& f4) const{return Dw_->gamma5(f4);}
  
  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;
  const Field gamma5(const Field&) const;

  /*!
   * @brief Calculates the \f$L_+(m)\f$
   */
  const Field proj_p(const Field& f4) const{
    Field w4=f4;
    w4 += Dw_->gamma5(f4);
    w4 *=0.5;
    return w4;
  }
    
  /*!
   * @brief Calculates the \f$L_-(m)\f$
   */
  const Field proj_m(const Field& f4) const{
    Field w4=f4;
    w4 -= Dw_->gamma5(f4);
    w4 *=0.5;
    return w4;
  }

  const Field Bproj(const Field& v5d) const;
  const Field Bproj_dag(const Field& v4d) const;
  const Field R5(const Field&) const;
  const Field R5g5(const Field&) const;

  const Field md_force( const Field& eta,const Field& zeta) const;
};

namespace DomainWallFermions {

  std::vector<double> getOmega(int Ns, 
			       double lambda_min,
			       double lambda_max);

}

#endif
