/*!
 * @file dirac_DomainWallCore.hpp
 * @brief Declaration of class Dirac_DomainWallCore (5d operator)
 Time-stamp: <2013-12-05 11:22:59 noaki>
 */
#ifndef DIRAC_DOMAINWALLCORE_INCLUDED
#define DIRAC_DOMAINWALLCORE_INCLUDED

#include "dirac_wilson.hpp"
#include "domainWallUtils.hpp"

class Proj5D;

////////////////////////////////////////////////////////////////////////////
/*! @brief Defines the 5d Domain Wall operator
 *
 * Defines \f[D_{\rm dwf}(m_q) = \omega D_W(-m_0)(1+cL(m_q))+(1-L(m_q))\f]
 * with: \f[L(m) = P_+ L_+(m) + P_- L_-(m)\f]
 * using proj_p() e proj_m() methods
 *
 * Kernel is given by \f[H = \frac{\gamma_5 b D_W}{ 2 + c D_W } \f]
 */
class DomainWallCore {
private:
  const DomainWallParams& prms_;
  const DiracWilsonLike* Dw_; /*!< @brief Dirac Kernel - any WilsonLike op */ 
  const Proj5D& projP_;
  const Proj5D& projM_;

  int N5_,f4size_;

  const Field get4d(const Field& f5,int s)const;
  void set5d(Field& f5,const Field& f4,int s)const;
  void add5d(Field& f5,const Field& f4,int s)const;
  void add5d(Field& f5,const Field& f4_1,const Field& f4_2,int s)const;
  void mul5d(Field& f5,double fac,int s)const;

  void(DomainWallCore::*mult_core)(Field&,const Field&)const;
  void(DomainWallCore::*mult_dag_core)(Field&,const Field&)const;
public:
  DomainWallCore(const DomainWallParams& prms,const DiracWilsonLike* Dw,
		 const Proj5D& projP,const Proj5D& projM)
    :prms_(prms),Dw_(Dw),projP_(projP),projM_(projM),
     N5_(prms.N5_),f4size_(Dw->fsize()),
     mult_core(&DomainWallCore::mult_full),
     mult_dag_core(&DomainWallCore::mult_dag_full){assert(Dw_);}

  DomainWallCore(const DomainWallParams& prms,const DiracWilsonLike* Dw,
		 const Proj5D& projP,const Proj5D& projM,
		 DWF::EvenOdd_tag)
    :prms_(prms),Dw_(Dw),projP_(projP),projM_(projM),
     N5_(prms.N5_),f4size_(Dw->fsize()),
     mult_core(&DomainWallCore::mult_offdiag),
     mult_dag_core(&DomainWallCore::mult_dag_offdiag){assert(Dw_);}
  
  /*@! brief mult in the heavy M0 limit */
  const Field mult_hop5(const Field& f5)const;    
  const Field mult_hop5_inv(const Field& f5)const;
  const Field mult_hop5_dag(const Field& f5)const;
  const Field mult_hop5_dinv(const Field& f5)const;

  //  void Dminus(Field&,const Field&)const;
  const Field mult(const Field& f5)const;
  const Field mult_dag(const Field& f5)const;
  void mult_offdiag(Field&,const Field&)const;    /*! @brief  -kpp*D*f    */
  void mult_full(Field&,const Field&)const;       /*! @brief  (1-kpp*D)*f */
  void mult_dag_offdiag(Field&,const Field&)const;/*! @brief -kpp*D^dag*f    */
  void mult_dag_full(Field&,const Field&)const;   /*! @brief (1-kpp*D^dag)*f */

  void md_force_p(Field&,const Field&,const Field&)const;
  void md_force_m(Field&,const Field&,const Field&)const;
  const Field md_force(const Field& phi,const Field& psi)const;
};

#endif
