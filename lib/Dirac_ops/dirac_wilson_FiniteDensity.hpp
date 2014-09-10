/*! @file dirac_wilson_FiniteDensity.hpp
 * @brief Dirac_Wilson_FiniteDensity class, which accommodates 
    finite chemical potential (the normal site-indexing version)
 Time-stamp: <2014-08-27 21:57:44 noaki>
 */
#ifndef DIRAC_WILSON_FINITEDENSITY_INCLUDED
#define DIRAC_WILSON_FINITEDENSITY_INCLUDED

#include "dirac_WilsonLikeFiniteDensity.hpp"
#include "dirac_wilson.hpp"


class Dirac_Wilson_FiniteDensity: public DiracWilsonLikeFiniteDensity{
private:
  const Field* const u_;
  Dirac_Wilson Dw_;
  double fgcty_; 

  std::valarray<size_t> sTmax_;
  std::valarray<size_t> sTmin_;
  void init_Tbdry();
public:
  Dirac_Wilson_FiniteDensity(double mass,double mu,const Field* u)
    :u_(u),Dw_(mass,u),fgcty_(exp(mu*CommonPrms::Lt())){init_Tbdry();}

  Dirac_Wilson_FiniteDensity(const XML::node& node,const Field* u)
    :u_(u),Dw_(node,u){
    double mu;
    XML::read(node, "mu",mu,MANDATORY);
    fgcty_= exp(mu*CommonPrms::Lt());

    init_Tbdry();
  }

  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;
  const Field md_force(const Field&,const Field&)const;

  const Field gamma5(const Field& f) const{ return Dw_.gamma5(f);}
  const Field* getGaugeField_ptr()const{ return u_;}  

  const Field mult_Ds(const Field&)const;
  const Field mult_Dtp(const Field&)const;
  const Field mult_Dtm(const Field&)const;

  size_t fsize() const{ return Dw_.fsize();}
  size_t gsize() const{ return Dw_.gsize();}
  double getMass()const{ return Dw_.getMass();}

#ifdef IBM_BGQ_WILSON
  void mult_ptr(double*,double* const )const;
  void mult_dag_ptr(double*,double* const )const;  
#endif
};

#endif
