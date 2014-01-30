/*! @file dirac_BFM_DomainWall_4D_eo.cpp
 *  @brief Methods of Dirac_BFM_DomainWall_4D_eo class

 * Time-stamp: <2014-01-30 14:13:20 neo>
 */
#include "dirac_BFM_DomainWall_4D_eo.hpp"
#include "Fields/field_expressions.hpp"

using namespace std;
using namespace FieldExpression;

#include "include/timings.hpp"
#include "include/messages_macros.hpp"
  
const Field Dirac_BFM_DomainWall_4D_eo::mult(const Field& f)const{
  // Original
  /*
  long double timing;
  FINE_TIMING_START(timing);

  Field sol5(fsize_*N5_);


  invDpv_->invert(sol5,invD_->mult(Bproj_dag(f)));// Dpv_^-1 
  FINE_TIMING_END(timing);
  _Message(TIMING_VERB_LEVEL,"[Timing] 4d DWF_BFM mult :"<<timing<<"\n"); 
  return Bproj(sol5);
  */

  ///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // WARNING: UNTESTED but it is similar to the mult_inv routine, just PV<->Op change
  ///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  long double timing, time_preparation, time_unprec;
  FINE_TIMING_START(timing);
  FINE_TIMING_START(time_preparation);
  //const ffmt_t ff_nex_(Nvol_,N5_);
  //Field sol5(fsize_*N5_);
  FermionField FField(Nvol_*N5_);
 
  // Using the fact the BFM assumes SiteIndexEO
  // so the vector is stored as
  // Vec[ s=0 Even | s=0 Odd | s=1 Even | s=1 Odd | s=2 Even | s=2 Odd | ... ]
  // where s is in the 5th dimension
  FField.data = Bproj_dag(f);
  /*
  Fermion_t* fm1 = BFM_Op_->LoadFullSource(FField);
  FINE_TIMING_END(time_preparation);
  FINE_TIMING_START(time_unprec);
  Fermion_t* fm2 = BFM_Op_->mult_unprec_base(fm1);
  FINE_TIMING_END(time_unprec);
  fm1 = BFM_Op_PV_->mult_inv_4d_base(fm2);
  BFM_Op_PV_->GetFullSolution(FField);
  FINE_TIMING_END(timing);
  */
  _Message(TIMING_VERB_LEVEL,"[Timing] 4d DWF_BFM mult_inv :"<<timing<<"\n"); 
  _Message(TIMING_VERB_LEVEL,"[Timing] 4d DWF_BFM mult_inv preparation :"<<time_preparation<<"\n"); 
  _Message(TIMING_VERB_LEVEL,"[Timing] 4d DWF_BFM mult_inv unprec :"<<time_unprec<<"\n"); 
  return Bproj(FField.data);
}

const Field Dirac_BFM_DomainWall_4D_eo::mult_inv(const Field& f)const{
  long double timing, time_preparation, time_unprec;
  FINE_TIMING_START(timing);
  FINE_TIMING_START(time_preparation);

  FermionField FField(Nvol_*N5_);
 
  // Using the fact the BFM assumes SiteIndexEO
  // so the vector is stored as
  // Vec[ s=0 Even | s=0 Odd | s=1 Even | s=1 Odd | s=2 Even | s=2 Odd | ... ]
  // where s is in the 5th dimension
  // which is the way of storing vectors in BFM
  FField.data = Bproj_dag(f);

  Fermion_t* fm1 = BFM_Op_PV_->LoadFullSource(FField);
  FINE_TIMING_END(time_preparation);
  FINE_TIMING_START(time_unprec);
  Fermion_t* fm2 = BFM_Op_PV_->mult_unprec_base(fm1);
  FINE_TIMING_END(time_unprec);
  fm1 = BFM_Op_->mult_inv_4d_base(fm2);
  BFM_Op_->GetFullSolution(FField);
  FINE_TIMING_END(timing);

  
  _Message(TIMING_VERB_LEVEL,"[Timing] 4d DWF_BFM mult_inv :"<<timing<<"\n"); 
  _Message(TIMING_VERB_LEVEL,"[Timing] 4d DWF_BFM mult_inv preparation :"<<time_preparation<<"\n"); 
  _Message(TIMING_VERB_LEVEL,"[Timing] 4d DWF_BFM mult_inv unprec :"<<time_unprec<<"\n"); 
  return Bproj(FField.data);


}

const Field Dirac_BFM_DomainWall_4D_eo::mult_dag(const Field& f)const{ 
  return gamma5(mult(gamma5(f)));}

const Field Dirac_BFM_DomainWall_4D_eo::mult_dag_inv(const Field& f)const{ 
  return gamma5(mult_inv(gamma5(f)));}


const Field Dirac_BFM_DomainWall_4D_eo::Bproj(const Field& f5) const{ 
  Field f4(fsize_);
#pragma omp parallel for
  for(int site=0; site<Nvol_; ++site)
    d5_.BprojCore(f4.getaddr(ff_.index(0,site)),
		  const_cast<Field&>(f5).getaddr(ff_.index(0,site,0)),
		  const_cast<Field&>(f5).getaddr(ff_.index(0,site,N5_-1)));
  return f4;
}

const Field Dirac_BFM_DomainWall_4D_eo::Bproj_dag(const Field& f4) const{
  Field f5(fsize_*N5_);
#pragma omp parallel for
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


Dirac_BFM_DomainWall_4D_eo::~Dirac_BFM_DomainWall_4D_eo(){
  CCIO::cout << "Destroying Dirac_BFM_DomainWall_4D_eo\n";
}
