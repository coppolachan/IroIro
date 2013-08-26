/*!--------------------------------------------------------------------------
 *@file dirac_DomainWall_EvenOdd.cpp
 *
 *@brief Definition of class methods for Dirac_optimalDomainWall_EvenOdd (5d op)

 Time-stamp: <2013-08-22 16:26:24 cossu>
 *-------------------------------------------------------------------------*/
#include "dirac_DomainWall_EvenOdd.hpp"
#include "Communicator/comm_io.hpp"
#include<stdlib.h>
#include<stdio.h>
#include<cassert>
#include<math.h>

#include "include/timings.hpp"
#include "include/messages_macros.hpp"

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
  Field w(Deo_.fsize()); // just slightly faster (but only BGQ)
  //Doe_->mult_hop(w,f);
  
  double* f_ptr = const_cast<Field&>(f).getaddr(0);
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
  timeval start_, end_;
  gettimeofday(&start_,NULL);

  Field w(Deo_.fsize());
  double* f_ptr = const_cast<Field&>(f).getaddr(0);
#pragma omp parallel
  {
    Deo_.mult_hop_dag_omp(w,f_ptr);
  }
  gettimeofday(&end_,NULL);
  multdag_timer += (end_.tv_sec -start_.tv_sec)*1000.0;
  multdag_timer += (end_.tv_usec -start_.tv_usec)/1000.0;   // us to ms
#else
  Field w(f);
  w -= mult_oe_dag(mult_eo_dag(f));
#endif   
  return w;
}

void Dirac_optimalDomainWall_EvenOdd::
md_force_eo(Field& fce, const Field& eta,const Field& zeta) const{
  Deo_.md_force_p(fce,eta,mult_ee_dinv(zeta));
  Doe_.md_force_m(fce,eta,mult_ee_dinv(zeta));
}

void Dirac_optimalDomainWall_EvenOdd::
md_force_oe(Field& fce, const Field& eta,const Field& zeta) const{
  Doe_.md_force_p(fce,eta,mult_oo_dinv(zeta));
  Deo_.md_force_m(fce,eta,mult_oo_dinv(zeta));
}

const Field Dirac_optimalDomainWall_EvenOdd::
md_force(const Field& eta,const Field& zeta) const{
  Field fce(Deo_.gsize());
  long double timing;
  FINE_TIMING_START(timing);

  md_force_eo(fce,mult_oe(eta),zeta);

  FINE_TIMING_END(timing);
  _Message(DEBUG_VERB_LEVEL, "[Timing] - Dirac_optimalDomainWall_EvenOdd::md_force"
           << " - md_force_eo timing = "
           << timing << std::endl);

  FINE_TIMING_START(timing);
  md_force_oe(fce,eta,mult_eo_dag(zeta));
  FINE_TIMING_END(timing);
  _Message(DEBUG_VERB_LEVEL, "[Timing] - Dirac_optimalDomainWall_EvenOdd::md_force"
           << " - md_force_oe timing = "
           << timing << std::endl);

  FINE_TIMING_START(timing);
  fce *= 0.5;
  FINE_TIMING_END(timing);
  _Message(DEBUG_VERB_LEVEL, "[Timing] - Dirac_optimalDomainWall_EvenOdd::md_force"
           << " - 0.5 mult timing = "
           << timing << std::endl);
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
void Dirac_optimalDomainWall_EvenOdd::solve_ms_eo(vector<Field>& xq,
						  const Field& b,
						  SolverOutput& SO, 
						  const vector<double>& sigma,
						  int MaxIter, 
						  double target_prec) const{
  Deo_.solve_ms_eo_5d(xq,b,SO,sigma,MaxIter,target_prec);
}

#endif
///////////////////////////////////////////////////////////////////////


