/*! @file dirac_DomainWall_4D_eoSolv.cpp
 *  @brief Methods of Dirac_optimalDomainWall_4D_eoSolv class
 */
#include "dirac_DomainWall_4D_eoSolv.hpp"
#include "Fields/field_expressions.hpp"

using namespace std;
using namespace FieldExpression;

const Field Dirac_optimalDomainWall_4D_eoSolv::mult(const Field& f)const{
  Field sol5(fsize_*N5_);
  invDpv_->invert(sol5,invD_->mult(Bproj_dag(f)));// Dpv_^-1 
  return Bproj(sol5);
}

const Field Dirac_optimalDomainWall_4D_eoSolv::mult_inv(const Field& f)const{
  Field sol5(fsize_*N5_);
  invD_->invert(sol5,invDpv_->mult(Bproj_dag(f)));// Dodw_^-1
  return Bproj(sol5);
}

const Field Dirac_optimalDomainWall_4D_eoSolv::mult_dag(const Field& f)const{ 
  return gamma5(mult(gamma5(f)));}

const Field Dirac_optimalDomainWall_4D_eoSolv::mult_dag_inv(const Field& f)const{ 
  return gamma5(mult_inv(gamma5(f)));}


const Field Dirac_optimalDomainWall_4D_eoSolv::Bproj(const Field& f5) const{ 
  Field f4(fsize_);
  for(int site=0; site<Nvol_; ++site)
    d5_.BprojCore(f4.getaddr(ff_.index(0,site)),
		  const_cast<Field&>(f5).getaddr(ff_.index(0,site,0)),
		  const_cast<Field&>(f5).getaddr(ff_.index(0,site,N5_-1)));
  return f4;
}

const Field Dirac_optimalDomainWall_4D_eoSolv::Bproj_dag(const Field& f4) const{
  Field f5(fsize_*N5_);
  for(int site=0; site<Nvol_; ++site)
    d5_.BprojCore_dag(f5.getaddr(ff_.index(0,site,0)),
		      f5.getaddr(ff_.index(0,site,N5_-1)),
		      const_cast<Field&>(f4).getaddr(ff_.index(0,site)));
  return f5;
}

const Field Dirac_optimalDomainWall_4D_eoSolv::gamma5(const Field& f)const{ 
  Field w(fsize_);
  for(int site=0; site<Nvol_; ++site)
    dm_.gamma5core(w.getaddr(ff_.index(0,site)),
		   const_cast<Field&>(f).getaddr(ff_.index(0,site)));
  return w;
}
/*
const Field Dirac_optimalDomainWall_4D_eoSolv::gamma5(const Field& f)const{ 
  Field w(D_->f4size());
  D_->gamma5_4d(w,f);
  return w;
}
*/
