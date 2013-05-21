/*!--------------------------------------------------------------------------
 *@file dirac_DomainWall_EvenOdd.cpp
 *
 *@brief Definition of class methods for Dirac_optimalDomainWall_EvenOdd (5d op)
 Time-stamp: <2013-05-21 09:34:21 noaki>
 *-------------------------------------------------------------------------*/
#include "dirac_DomainWall_EvenOdd.hpp"
#include "Communicator/comm_io.hpp"
#include<stdlib.h>
#include<stdio.h>
#include<cassert>
#include<math.h>

using namespace std;

const Field Dirac_optimalDomainWall_EvenOdd::mult_ee(const Field& f)const{
  return Deo_.mult_hop5(f);
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_oo(const Field& f)const{
  return Deo_.mult_hop5(f);
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_ee_inv(const Field& f)const{
  return Deo_.mult_hop5_inv(f);
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_ee_dinv(const Field& f)const{
  return Deo_.mult_hop5_dinv(f);
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_oo_inv(const Field& f)const{
  return Deo_.mult_hop5_inv(f);
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_oo_dinv(const Field& f)const{
  return Deo_.mult_hop5_dinv(f);
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_eo(const Field& f)const{
  return Deo_.mult_hop5_inv(Deo_.mult(f));
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_oe(const Field& f)const{
  return Doe_.mult_hop5_inv(Doe_.mult(f));
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_eo_dag(const Field& f)const{
  return Doe_.mult_dag(Deo_.mult_hop5_dinv(f));
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_oe_dag(const Field& f)const{
  return Deo_.mult_dag(Doe_.mult_hop5_dinv(f));
}

const Field Dirac_optimalDomainWall_EvenOdd::mult(const Field& f) const{
#ifdef  IBM_BGQ_WILSON
  double* f_ptr = const_cast<Field&>(f).getaddr(0);
  Field w(Deo_.fsize()); // just slightly faster (but only BGQ)

#pragma omp parallel
  {
    Deo_.mult_hop_omp(w,f_ptr);
  }
#else
  Field w(f);
  w -= mult_eo(mult_oe(f));
#endif
  return w;
}

const Field Dirac_optimalDomainWall_EvenOdd::mult_dag(const Field& f) const{
#ifdef IBM_BGQ_WILSON  
  double* f_ptr = const_cast<Field&>(f).getaddr(0);
  Field w(Deo_.fsize());

#pragma omp parallel
  {
    Deo_.mult_hop_dag_omp(w,f_ptr);
  }
#else
  Field w(f);
  w -= mult_oe_dag(mult_eo_dag(f));
#endif   
  return w;
}

void Dirac_optimalDomainWall_EvenOdd::
md_force_eo(Field& fce, const Field& eta,const Field& zeta) const{
#ifdef IBM_BGQ_WILSON  
  Deo_.md_force_p_BGQ(fce,eta,mult_ee_dinv(zeta));
  Doe_.md_force_m_BGQ(fce,eta,mult_ee_dinv(zeta));
#else
  Deo_.md_force_p(fce,eta,mult_ee_dinv(zeta));
  Doe_.md_force_m(fce,eta,mult_ee_dinv(zeta));
#endif
}

void Dirac_optimalDomainWall_EvenOdd::
md_force_oe(Field& fce, const Field& eta,const Field& zeta) const{
#ifdef IBM_BGQ_WILSON  
  Doe_.md_force_p_BGQ(fce,eta,mult_oo_dinv(zeta));
  Deo_.md_force_m_BGQ(fce,eta,mult_oo_dinv(zeta));
#else
  Doe_.md_force_p(fce,eta,mult_oo_dinv(zeta));
  Deo_.md_force_m(fce,eta,mult_oo_dinv(zeta));
#endif
}

const Field Dirac_optimalDomainWall_EvenOdd::
md_force(const Field& eta,const Field& zeta) const{
  Field fce(Deo_.gsize());
  md_force_eo(fce,mult_oe(eta),zeta);
  md_force_oe(fce,eta,mult_eo_dag(zeta));
  fce *= 0.5;
  return fce;
}

///////////////////////////////////////////////////////////////////////
#ifdef IBM_BGQ_WILSON
// OPTIMIZED LIBRARIES
void  Dirac_optimalDomainWall_EvenOdd::solve_eo(Field& out, 
						const Field& in, 
						SolverOutput& SO, 
						int Niter, 
						double stop_cond) const{
  Deo_.solve_eo_5d(out,in,SO,Niter,stop_cond);
}
void Dirac_optimalDomainWall_EvenOdd::solve_ms_eo(prop_t& xq,
						  const Field& b,
						  SolverOutput& SO, 
						  const vector<double>& sigma,
						  int MaxIter, 
						  double GoalPrecision) const{
  Deo_.solve_ms_eo_5d(xq,b,SO,sigma,MaxIter,GoalPrecision);
}

#endif
///////////////////////////////////////////////////////////////////////


