/*!
 * @file dirac_DomainWall.hpp
 *
 * @brief Declaration of class Dirac_optimalDomainWall (5d operator)
 */
#ifndef DIRAC_OPTIMALDOMAINWALL_INCLUDED
#define DIRAC_OPTIMALDOMAINWALL_INCLUDED

#include <string.h>
#include "include/field.h"
#include "include/pugi_interface.h"
#include "Dirac_ops/dirac_wilson.h"

namespace DomainWallFermions {

  std::vector<double> getOmega(int Ns, 
			       double lambda_min,
			       double lambda_max);

}

enum {Standard, PauliVillars};

/*!
 * @brief Container for parameter of the 5d Optimal Domain Wall operator
 *
 * Parameters for \f[D_{\rm dwf}(m_q) = \omega D_W(-m_0)(1+cL(m_q))+(1-L(m_q))\f]
 * 
 *
 */
struct Dirac_optimalDomainWall_params{
  double N5dim_;
  double b_, c_;
  double mq_;/*!< @brief Bare quark mass \f$m_q\f$ */
  std::vector<double> omega_;/*!< @brief Diagonal matrix
			      * \f$\omega = \{\omega_s\}\f$
			      */
  std::vector<double> bs_, cs_;
  std::vector<double> dp_, dm_;

  Dirac_optimalDomainWall_params(XML::node DWF_node){
    XML::read(DWF_node, "N5d", N5dim_, MANDATORY);
    XML::read(DWF_node, "b", b_, MANDATORY);
    XML::read(DWF_node, "c", c_, MANDATORY);
    XML::read(DWF_node, "mass", mq_, MANDATORY);
    omega_.resize(N5dim_);
    XML::node ApproxNode = DWF_node.child("approximation");
    if (ApproxNode !=NULL) {
      const char* Approx_name = ApproxNode.attribute("name").value();
      if (!strcmp(Approx_name, "Zolotarev")){
	double lambda_min, lambda_max;
	XML::read(ApproxNode, "lambda_min", lambda_min); 
	XML::read(ApproxNode, "lambda_max", lambda_max); 
	omega_ = DomainWallFermions::getOmega(N5dim_,lambda_min,lambda_max);
      }
      if (!strcmp(Approx_name, "Tanh")){
	for (int s = 0; s < N5dim_; ++s) omega_[s] = 1.0;
      }
    } else {
      abort();
    }

    bs_.resize(N5dim_);
    cs_.resize(N5dim_);
    dp_.resize(N5dim_);
    dm_.resize(N5dim_);

  }

  Dirac_optimalDomainWall_params(const double b,
				 const double c,
				 const double mass,
				 const std::vector<double> omega){
    b_ = b;
    c_ = c;
    mq_ = mass;
    omega_ = omega;
    bs_.resize(omega.size());
    cs_.resize(omega.size());
    dp_.resize(omega.size());
    dm_.resize(omega.size());

    // this line must be replaced
    double M0_=-1.6;  
    //
    for (int s = 0; s < omega.size(); ++s) {
      bs_[s] = (b*omega[s]+c)/2.0;
      cs_[s] = (b*omega[s]-c)/2.0;
      dp_[s] = bs_[s]*(4.0+M0_)+1.0;
      dm_[s] = 1.0-cs_[s]*(4.0+M0_);
      //      std::cout << s << " " << dp_[s] << " " << dm_[s] << std::endl;
    }
  }

  Dirac_optimalDomainWall_params(const Dirac_optimalDomainWall_params& Par,
				 int Type = 0){
    switch (Type){
    case Standard:
      b_ = Par.b_;
      c_ = Par.c_;
      mq_ = Par.mq_;
      omega_ = Par.omega_;
      bs_ = Par.bs_;
      cs_ = Par.cs_;
      dp_ = Par.dp_;
      dm_ = Par.dm_;
      break;
    case PauliVillars:
      b_ = Par.b_;
      c_ = Par.c_;
      mq_ = 1.0;
      omega_ = Par.omega_;
      bs_ = Par.bs_;
      cs_ = Par.cs_;
      dp_ = Par.dp_;
      dm_ = Par.dm_;
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
  Dirac_optimalDomainWall_params Params;

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
    
    for (int s = 0; s < N5_; ++s) {
      Params.bs_[s] = (Params.b_*Params.omega_[s]+Params.c_)/2.0;
      Params.cs_[s] = (Params.b_*Params.omega_[s]-Params.c_)/2.0;
      Params.dp_[s] = Params.bs_[s]*(4.0+M0_)+1.0;
      Params.dm_[s] = 1.0-Params.cs_[s]*(4.0+M0_);
      //      std::cout << s << " " << dp_[s] << " " << dm_[s] << std::endl;
    }
    
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

 
  Dirac_optimalDomainWall(const double b,
			  const double c,
			  const double mq,
			  const std::vector<double>& omega,
			  const Dirac_Wilson* Kernel)
    :Params(b,c,mq,omega),
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
    Field w4 = Dw_->proj_p(f4);
    return w4;
  }
  //  const Field proj_p(const Field& f4) const{
  //    Field w4=f4;
  //    w4 += Dw_->gamma5(f4);
  //    w4 *=0.5;
  //    return w4;
  //  }
    
  /*!
   * @brief Calculates the \f$L_-(m)\f$
   */
  const Field proj_m(const Field& f4) const{
    Field w4 = Dw_->proj_m(f4);
    return w4;
  }
  //  const Field proj_m(const Field& f4) const{
  //    Field w4=f4;
  //    w4 -= Dw_->gamma5(f4);
  //    w4 *=0.5;
  //    return w4;
  //  }

  const Field Bproj(const Field& v5d) const;
  const Field Bproj_dag(const Field& v4d) const;
  const Field R5(const Field&) const;
  const Field R5g5(const Field&) const;

  const Field md_force( const Field& eta,const Field& zeta) const;
};


#endif
