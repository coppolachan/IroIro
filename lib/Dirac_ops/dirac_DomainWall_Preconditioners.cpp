/*! 
  @file dirac_DomainWall_Preconditioners.cpp
  
  @brief Definition of the Preconditioners for DomainWall operator

*/
#include "dirac_DomainWall.hpp"
#include "Fields/field_expressions.hpp"
#include "include/field.h"

//---------------------------------------------------------------------------
const Field Dirac_optimalDomainWall::NoPrecond::mult(const Field& f5) const{
  return DWF_->mult(f5);
}

const Field Dirac_optimalDomainWall::NoPrecond::mult_dag(const Field& f5) const{
  return DWF_->mult_dag(f5);
}

const Field Dirac_optimalDomainWall::NoPrecond::left(const Field& f5) const{
  return f5;
}
const Field Dirac_optimalDomainWall::NoPrecond::right(const Field& f5) const{
  return f5;
}
const Field Dirac_optimalDomainWall::NoPrecond::left_dag(const Field& f5) const{
  return f5;
}
const Field Dirac_optimalDomainWall::NoPrecond::right_dag(const Field& f5) const{
  return f5;
}

//-----------------------------------------------------------------------------
const Field Dirac_optimalDomainWall::LUPrecond::mult(const Field& f5) const{
  using namespace FieldExpression;
  assert(f5.size()==DWF_->fsize_);
  Field w5(DWF_->fsize_);
  w5 = DWF_->mult(f5);

  //  Field t5(DWF_->fsize_);
  //  t5 = LU_inv(w5);
  //  Field u5(DWF_->fsize_);
  //  u5 = LU(t5);
  //  CCIO::cout << "LU before " << w5.norm() << std::endl;
  //  CCIO::cout << "LUinv     " << t5.norm() << std::endl;
  //  CCIO::cout << "LU LUinv  " << u5.norm() << std::endl;

  return LU_inv(w5);
}

const Field Dirac_optimalDomainWall::LUPrecond::mult_dag(const Field& f5) const{
  assert(f5.size()==DWF_->fsize_);
  Field t5(DWF_->fsize_);
  t5 = LU_dag_inv(f5);

  //  Field u5(DWF_->fsize_);
  //  u5 = LU_dag(t5);
  //  CCIO::cout << "LUdag before    " << f5.norm() << std::endl;
  //  CCIO::cout << "LUdaginv        " << t5.norm() << std::endl;
  //  CCIO::cout << "LUdag LUdaginv  " << u5.norm() << std::endl;

  return DWF_->mult_dag(t5); 
}

const Field Dirac_optimalDomainWall::LUPrecond::LU(const Field& f5) const{
 using namespace FieldExpression;
  
  assert(f5.size()==DWF_->fsize_);
  Field w5(DWF_->fsize_);

  for (int s=0; s<DWF_->N5_; ++s) {
    Field v = DWF_->get4d(f5,s);
    v *= DWF_->Params.dp_[s];
    DWF_->set5d(w5,v,s);
  }
  for (int s=0; s<DWF_->N5_-1; ++s) {
    Field v = DWF_->proj_m(DWF_->get4d(f5,s+1));
    v *= -DWF_->Params.dm_[s];
    DWF_->add5d(w5,v,s);
  }
  for (int s=1; s<DWF_->N5_; ++s) {
    Field v = DWF_->proj_p(DWF_->get4d(f5,s-1));
    v *= -DWF_->Params.dm_[s];
    DWF_->add5d(w5,v,s);
  }

  Field v = DWF_->proj_p(DWF_->get4d(f5,DWF_->N5_-1));
  v *= DWF_->Params.mq_*DWF_->Params.dm_[0];
  DWF_->add5d(w5,v,0);

  v = DWF_->proj_m(DWF_->get4d(f5,0));
  v *= DWF_->Params.mq_*DWF_->Params.dm_[DWF_->N5_-1];
  DWF_->add5d(w5,v,DWF_->N5_-1);

  return w5;
}

const Field Dirac_optimalDomainWall::LUPrecond::LU_dag(const Field& f5) const{
 using namespace FieldExpression;
  
  assert(f5.size()==DWF_->fsize_);
  Field w5(DWF_->fsize_);

  for (int s=0; s<DWF_->N5_; ++s) {
    Field v = DWF_->get4d(f5,s);
    v *= DWF_->Params.dp_[s];
    DWF_->set5d(w5,v,s);
  }
  for (int s=0; s<DWF_->N5_-1; ++s) {
    Field v = DWF_->proj_p(DWF_->get4d(f5,s+1));
    v *= -DWF_->Params.dm_[s+1];
    DWF_->add5d(w5,v,s);
  }
  for (int s=1; s<DWF_->N5_; ++s) {
    Field v = DWF_->proj_m(DWF_->get4d(f5,s-1));
    v *= -DWF_->Params.dm_[s-1];
    DWF_->add5d(w5,v,s);
  }

  Field v = DWF_->proj_m(DWF_->get4d(f5,DWF_->N5_-1));
  v *= DWF_->Params.mq_*DWF_->Params.dm_[DWF_->N5_-1];
  DWF_->add5d(w5,v,0);

  v = DWF_->proj_p(DWF_->get4d(f5,0));
  v *= DWF_->Params.mq_*DWF_->Params.dm_[0];
  DWF_->add5d(w5,v,DWF_->N5_-1);

  return w5;
}

const Field Dirac_optimalDomainWall::LUPrecond::LU_inv(const Field& f5) const{
 using namespace FieldExpression;
  
  assert(f5.size()==DWF_->fsize_);
  //  Field w5(DWF_->fsize_);
  Field w5(f5);

  for (int s=1; s<DWF_->N5_; ++s) {
    Field lpf = DWF_->proj_p(DWF_->get4d(w5,s-1));
    lpf *= (DWF_->Params.dm_[s]/DWF_->Params.dp_[s-1]);
    DWF_->add5d(w5,lpf,s);
  }
  for (int s=0; s<DWF_->N5_-1; ++s) {
    Field ey = DWF_->proj_m(DWF_->get4d(w5,s));
    ey *= - DWF_->Params.mq_*DWF_->Params.es_[s];
    DWF_->add5d(w5,ey,DWF_->N5_-1);
  }

  Field v = DWF_->get4d(w5,DWF_->N5_-1);
  v *= 1.0 / (DWF_->Params.dp_[DWF_->N5_-1] +
	      DWF_->Params.mq_*
	      DWF_->Params.dm_[DWF_->N5_-2]*DWF_->Params.es_[DWF_->N5_-2]);
  DWF_->set5d(w5,v,DWF_->N5_-1);
  for (int s=DWF_->N5_-2; s>=0; --s) {
    Field lmf = DWF_->proj_m(DWF_->get4d(w5,s+1));
    lmf *= DWF_->Params.dm_[s];
    DWF_->add5d(w5,lmf,s);
    Field fy = DWF_->proj_p(DWF_->get4d(w5,DWF_->N5_-1));
    fy *= - DWF_->Params.mq_*DWF_->Params.fs_[s];
    DWF_->add5d(w5,fy,s);
    v = DWF_->get4d(w5,s);
    v *= 1.0/DWF_->Params.dp_[s];
    DWF_->set5d(w5,v,s);
  }
  return w5;

}

const Field Dirac_optimalDomainWall::LUPrecond::LU_dag_inv(const Field& f5) const{
   
 assert(f5.size()==DWF_->fsize_);
  Field t5(DWF_->fsize_);
  t5 = f5;
  
  // LU preconditioning : ((LU)^T)^-1 = (U^T L^T)^-1 = (L^T)^-1 (U^T)^-1
  Field v = DWF_->get4d(t5,0);
  v *= 1.0/DWF_->Params.dp_[0];
  DWF_->set5d(t5,v,0);

  for (int s=1; s<DWF_->N5_-1; ++s) {
    Field lmf = DWF_->proj_m(DWF_->get4d(t5,s-1));
    lmf *= DWF_->Params.dm_[s-1];
    DWF_->add5d(t5,lmf,s);
    v = DWF_->get4d(t5,s);
    v *= 1.0/DWF_->Params.dp_[s];
    DWF_->set5d(t5,v,s);
  }

  v = DWF_->proj_m(DWF_->get4d(t5,DWF_->N5_-2));
  v *= DWF_->Params.dm_[DWF_->N5_-2];
  DWF_->add5d(t5,v,DWF_->N5_-1);
  for (int s=0; s<DWF_->N5_-1; ++s) {
    Field fy = DWF_->proj_p(DWF_->get4d(t5,s));
    fy *= -DWF_->Params.mq_*DWF_->Params.fs_[s];
    DWF_->add5d(t5,fy,DWF_->N5_-1);
  }
  v = DWF_->get4d(t5,DWF_->N5_-1);
  v *= 1.0 / (DWF_->Params.dp_[DWF_->N5_-1] +
	      DWF_->Params.mq_*
	      DWF_->Params.dm_[DWF_->N5_-2]*DWF_->Params.es_[DWF_->N5_-2]);
  DWF_->set5d(t5,v,DWF_->N5_-1);

  for (int s=DWF_->N5_-2; s>=0; --s) {
    Field lpf = DWF_->proj_p(DWF_->get4d(t5,s+1));
    lpf *= (DWF_->Params.dm_[s+1]/DWF_->Params.dp_[s]);
    DWF_->add5d(t5,lpf,s);
    Field ey = DWF_->proj_m(DWF_->get4d(t5,DWF_->N5_-1));
    ey *= - DWF_->Params.mq_*DWF_->Params.es_[s];
    DWF_->add5d(t5,ey,s);
  }

  return t5; 
}

const Field Dirac_optimalDomainWall::LUPrecond::left(const Field& f5) const{
  return LU_inv(f5);
  //  return f5;
}
const Field Dirac_optimalDomainWall::LUPrecond::right(const Field& f5) const{
  return f5;
  //  return LU_inv(f5);
}
const Field Dirac_optimalDomainWall::LUPrecond::left_dag(const Field& f5) const{
  return LU_dag_inv(f5);
  //  return f5;
}
const Field Dirac_optimalDomainWall::LUPrecond::right_dag(const Field& f5) const{
  return f5;
  //  return LU_dag_inv(f5);
}
//-------------------------------------------------------------------------
