/*!
 * @file dirac_DomainWall_4D_fullSolv.cpp
 * @brief Declaration of Dirac_optimalDomainWall_4D_fullSolv class 
 Time-stamp: <2013-05-31 07:15:18 noaki>
 */
#include "dirac_DomainWall_4D_fullSolv.hpp"
#include "Fields/field_expressions.hpp"

using namespace std;
using namespace FieldExpression;

void Dirac_optimalDomainWall_4D_fullSolv::mult_normal(Field& w,const Field& f)const{
  // D_dw(m)
  w = Dodw_->mult(Dodw_->Bproj_dag(f));
  // Dpv_^-1
  Field src = Dpv_->mult_dag_prec(Dpv_->left_prec(w));
  Field sol5(fsize_);
  SolverOutput monitor = slv_pv_->solve(sol5,src);
#if VERBOSITY > 0
  monitor.print();
#endif
  w = Dodw_->Bproj(Dpv_->right_prec(sol5));
}

// mult with Exact low-modes
void Dirac_optimalDomainWall_4D_fullSolv::mult_lmp(Field& w,const Field& f)const{
  Field f_h = lmh_->proj_high(f);
  w = (0.5*(1.0+mq_))*(f-f_h);
  w += (0.5*(1.0-mq_))*gamma5(lmh_->proj_appliedLow(f));
  Field v(fsize_);  
  mult_normal(v,f_h); 
  w += v;
}

void Dirac_optimalDomainWall_4D_fullSolv::mult_inv_normal(Field& w,const Field& f)const{
  // D_pv
  w = Dpv_->mult(Dodw_->Bproj_dag(f));
  // D_dw^-1
  Field src = Dodw_->mult_dag_prec(Dodw_->left_prec(w));
  Field sol5(Dodw_->fsize());
  SolverOutput monitor = slv_odw_->solve(sol5,src);
#if VERBOSITY > 0
  monitor.print();
#endif
 w = Dodw_->Bproj(Dodw_->right_prec(sol5));
}

// mult_inv with Exact low-modes
void Dirac_optimalDomainWall_4D_fullSolv::mult_inv_lmp(Field& w, const Field& f)const{ 
// This function is not tested, yet.
  Field f_h = lmh_->proj_high(f);  
  w = (0.5*(1.0+mq_)/mq_)*(f-f_h); 
  w += (0.5*(mq_-1.0)/mq_)*gamma5(lmh_->proj_appliedLow(f));
  Field v(fsize_);  
  mult_inv_normal(v,f_h); 
  w += v;
}

const Field Dirac_optimalDomainWall_4D_fullSolv::mult(const Field& f)const{
  Field w(fsize_);
  (this->*mult_core)(w,f);
  return w;
}

const Field Dirac_optimalDomainWall_4D_fullSolv::mult_dag(const Field& f)const{
  return gamma5(mult(gamma5(f)));}

const Field Dirac_optimalDomainWall_4D_fullSolv::mult_inv(const Field& f)const{
  Field w(fsize_);
  (this->*mult_inv_core)(w,f);
  return w;
}

const Field Dirac_optimalDomainWall_4D_fullSolv::mult_dag_inv(const Field& f)const{
  return gamma5(mult_dag(gamma5(f)));}

const Field Dirac_optimalDomainWall_4D_fullSolv::gamma5(const Field& f)const{ 
  Field w(Dodw_->f4size());
  Dodw_->gamma5_4d(w,f);
  return w;
}
