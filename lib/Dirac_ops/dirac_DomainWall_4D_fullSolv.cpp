/*!
 * @file dirac_DomainWall_4D_fullSolv.cpp
 * @brief Declaration of Dirac_optimalDomainWall_4D_fullSolv class 
 Time-stamp: <2013-04-23 17:54:12 noaki>
 */
#include "dirac_DomainWall_4D_fullSolv.hpp"
#include "Tools/randNum_MP.h"

using namespace std;

const Field Dirac_optimalDomainWall_4D_fullSolv::
mult(const Field& f)const{
  // D_dw(m)
  Field t5 = Dodw_->mult(Dodw_->Bproj_dag(f));
  // Dpv_^-1
  Field src = Dpv_->mult_dag_prec(Dpv_->left_prec(t5));
  Field sol5(Dodw_->fsize());
  SolverOutput monitor = slv_pv_->solve(sol5,src);
#if VERBOSITY > 0
  monitor.print();
#endif
  return Dodw_->Bproj(Dpv_->right_prec(sol5));
}

const Field Dirac_optimalDomainWall_4D_fullSolv::
mult_dag(const Field& f)const{
  // Dpv_^dag^-1
  Field src = Dpv_->right_dag_prec(Dodw_->Bproj_dag(f));
  Field sol5(Dodw_->fsize());
  SolverOutput monitor = slv_pv_->solve(sol5,src);
#if VERBOSITY > 0
  monitor.print();
#endif
  src = Dpv_->left_dag_prec(Dpv_->mult_prec(sol5));
  // D_dw^dag(m)
  return Dodw_->Bproj(Dodw_->mult_dag(src));
}

const Field Dirac_optimalDomainWall_4D_fullSolv::
mult_inv(const Field& f)const{
  // D_pv
  Field t5 = Dpv_->mult(Dodw_->Bproj_dag(f));
  // D_dw^-1
  Field src = Dodw_->mult_dag_prec(Dodw_->left_prec(t5));
  Field sol5(Dodw_->fsize());
  SolverOutput monitor = slv_odw_->solve(sol5,src);
#if VERBOSITY > 0
  monitor.print();
#endif
  CCIO::cout<<"sol5.norm()="<<(Dodw_->right_prec(sol5)).norm()<<std::endl;
  return Dodw_->Bproj(Dodw_->right_prec(sol5));
}

const Field Dirac_optimalDomainWall_4D_fullSolv::
mult_dag_inv(const Field& f)const{
  // D_dw^dag^-1
  Field src = Dodw_->right_dag_prec(Dodw_->Bproj_dag(f));
  Field sol5(Dodw_->fsize());
  SolverOutput monitor = slv_odw_->solve(sol5,src);
#if VERBOSITY > 0
  monitor.print();
#endif
  src = Dodw_->left_dag_prec(Dodw_->mult_prec(sol5));
  // D_pv^dag
  return Dodw_->Bproj(Dpv_->mult_dag(src));
}

const Field Dirac_optimalDomainWall_4D_fullSolv::gamma5(const Field& f)const{ 
  Field w(Dodw_->f4size());
  Dodw_->gamma5_4d(w,f);
  return w;
}
