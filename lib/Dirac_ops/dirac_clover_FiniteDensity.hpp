/*! @file dirac_wilson_FiniteDensity.hpp
 * @brief Dirac_Wilson_FiniteDensity class, which accommodates 
    finite chemical potential (the normal site-indexing version)
 Time-stamp: <2014-08-28 08:47:24 noaki>
 */
#ifndef DIRAC_CLOVER_FINITEDENSITY_INCLUDED
#define DIRAC_CLOVER_FINITEDENSITY_INCLUDED

/////// this is still in the beginning of working ///////

#include "dirac_wilson_FiniteDensity.hpp"

class Dirac_Clover_FiniteDensity: public DiracWilsonLikeFiniteDensity{
private:
  const Field* const u_;
  Dirac_Wilson_FiniteDensity Dwfd_;
public:
  Dirac_Clover_FiniteDensity(double mass,double mu,const Field* u)
    :Dwfd_(mass,mu,u),u_(u){}

  Dirac_Clover_FiniteDensity(const XML::node& node,const Field* u)
    :Dwfd_(node,u),u_(u){}

  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;
  const Field md_force(const Field&,const Field&)const;

  const Field gamma5(const Field& f) const{ return Dwfd_.gamma5(f);}
  const Field* getGaugeField_ptr()const{ return u_;}  

  const Field mult_Ds(const Field& f)const{ return Dwfd_.mult_Ds(f);}
  const Field mult_Dtp(const Field& f)const{ return Dwfd_.mult_Dtp(f);}
  const Field mult_Dtm(const Field& f)const{ return Dwfd_.mult_Dtm(f);}
  const Field mult_Ex(const Field& f)const{/* return  mult_sw(f); */}

  size_t fsize() const{ return Dwfd_.fsize();}
  size_t gsize() const{ return Dwfd_.gsize();}
  double getMass()const{ return Dwfd_.getMass();}

#ifdef IBM_BGQ_WILSON
  void mult_ptr(double*,double* const )const;
  void mult_dag_ptr(double*,double* const )const;  
#endif
};

#endif
