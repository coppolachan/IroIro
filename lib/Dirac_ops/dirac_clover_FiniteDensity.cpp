/*! @file dirac_clover_FiniteDensity.cpp
 * @brief memberfuncs of Dirac_Clover_FiniteDensity class
 Time-stamp: <2014-08-27 11:53:49 noaki>
*/
#include "dirac_clover_FiniteDensity.hpp"

using namespace std;

const Field Dirac_Clover_FiniteDensity::mult(const Field& f)const{
Field Df = Dwfd_.mult(f);
return Df;
}

const Field Dirac_Clover_FiniteDensity::mult_dag(const Field& f)const{
Field g5f = Dwfd_.mult_dag(f);

return g5f;
}

const Field Dirac_Clover_FiniteDensity::
md_force(const Field& eta,const Field& zeta)const{

  Field fce = Dwfd_.md_force(eta,zeta);
  return fce;
}
  
#ifdef IBM_BGQ_WILSON
void Dirac_Clover_FiniteDensity::
mult_ptr(double* w,double* const f)const{
  Dwfd_.mult_ptr(w,f);
}

void Dirac_Clover_FiniteDensity::
mult_dag_ptr(double* w,double* const f)const{
  Dwfd_.mult_dag_ptr(w,f);
}

#endif
