/*! @file dirac_BFM_DomainWall_4D_eo.cpp
 *  @brief Methods of Dirac_BFM_DomainWall_4D_eo class

 * Time-stamp: <2014-01-27 13:10:16 neo>
 */
#include "dirac_BFM_DomainWall_4D_eo.hpp"
#include "Fields/field_expressions.hpp"

using namespace std;
using namespace FieldExpression;

#include "include/timings.hpp"
#include "include/messages_macros.hpp"

const Field Dirac_BFM_DomainWall_4D_eo::mult(const Field& f)const{
  long double timing;
  FINE_TIMING_START(timing);

  Field sol5(fsize_*N5_);
  invDpv_->invert(sol5,invD_->mult(Bproj_dag(f)));// Dpv_^-1 
  FINE_TIMING_END(timing);
  _Message(TIMING_VERB_LEVEL,"[Timing] 4d DWF mult :"<<timing<<"\n"); 
  return Bproj(sol5);
}

const Field Dirac_BFM_DomainWall_4D_eo::mult_inv(const Field& f)const{
  long double timing;
  FINE_TIMING_START(timing);

  Field sol5(fsize_*N5_);
  invD_->invert(sol5,invDpv_->mult(Bproj_dag(f)));// Dodw_^-1
  FINE_TIMING_END(timing);
  _Message(TIMING_VERB_LEVEL,"[Timing] 4d DWF mult_inv :"<<timing<<"\n"); 
  return Bproj(sol5);


}

const Field Dirac_BFM_DomainWall_4D_eo::mult_dag(const Field& f)const{ 
  return gamma5(mult(gamma5(f)));}

const Field Dirac_BFM_DomainWall_4D_eo::mult_dag_inv(const Field& f)const{ 
  return gamma5(mult_inv(gamma5(f)));}


const Field Dirac_BFM_DomainWall_4D_eo::Bproj(const Field& f5) const{ 
  Field f4(fsize_);
  for(int site=0; site<Nvol_; ++site)
    d5_.BprojCore(f4.getaddr(ff_.index(0,site)),
		  const_cast<Field&>(f5).getaddr(ff_.index(0,site,0)),
		  const_cast<Field&>(f5).getaddr(ff_.index(0,site,N5_-1)));
  return f4;
}

const Field Dirac_BFM_DomainWall_4D_eo::Bproj_dag(const Field& f4) const{
  Field f5(fsize_*N5_);
  for(int site=0; site<Nvol_; ++site)
    d5_.BprojCore_dag(f5.getaddr(ff_.index(0,site,0)),
		      f5.getaddr(ff_.index(0,site,N5_-1)),
		      const_cast<Field&>(f4).getaddr(ff_.index(0,site)));
  return f5;
}

const Field Dirac_BFM_DomainWall_4D_eo::gamma5(const Field& f)const{ 
  Field w(fsize_);
  for(int site=0; site<Nvol_; ++site)
    dm_.gamma5core(w.getaddr(ff_.index(0,site)),
		   const_cast<Field&>(f).getaddr(ff_.index(0,site)));
  return w;
}
