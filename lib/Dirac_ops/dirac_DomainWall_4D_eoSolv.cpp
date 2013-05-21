/*! @file dirac_DomainWall_4D_eoSolv.cpp
 *  @brief Methods of Dirac_optimalDomainWall_4D_eoSolv class
 */
#include "dirac_DomainWall_4D_eoSolv.hpp"
#include "Fields/field_expressions.hpp"

using namespace std;

const Field Dirac_optimalDomainWall_4D_eoSolv::mult(const Field& f)const{
  Field sol5(D_->fsize());
  //invDpv_->invert(sol5,D_->mult(D_->Bproj_dag(f)));   // Dpv_^-1
  invDpv_->invert(sol5,invD_->mult(D_->Bproj_dag(f)));   // Dpv_^-1
  return D_->Bproj(sol5);
}

const Field Dirac_optimalDomainWall_4D_eoSolv::mult_dag(const Field& f)const{ 
  return gamma5(mult(gamma5(f)));}

const Field Dirac_optimalDomainWall_4D_eoSolv::mult_inv(const Field& f)const{
  Field sol5(D_->fsize());
  //invD_->invert(sol5,Dpv_->mult(D_->Bproj_dag(f)));   // Dodw_^-1
  invD_->invert(sol5,invDpv_->mult(D_->Bproj_dag(f)));   // Dodw_^-1
  return D_->Bproj(sol5);
}

const Field Dirac_optimalDomainWall_4D_eoSolv::mult_dag_inv(const Field& f)const{ 
  return gamma5(mult_inv(gamma5(f)));}

const Field Dirac_optimalDomainWall_4D_eoSolv::gamma5(const Field& f)const{ 
  Field w(D_->f4size());
  D_->gamma5_4d(w,f);
  return w;
}


