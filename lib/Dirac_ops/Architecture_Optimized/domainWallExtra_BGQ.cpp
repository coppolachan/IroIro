/*!--------------------------------------------------------------------------
 * @file domainWallExtra_BGQ.cpp
 *
 * @brief Definition of extra functions on BGQ for Dirac_DomainWall (5d op.)
 *Time-stamp: <2014-01-20 12:47:38 noaki>
 *-------------------------------------------------------------------------*/
#include "Tools/Architecture_Optimized/utils_BGQ.hpp"
#include "Dirac_ops/dirac_DomainWall.hpp"

void Dirac_DomainWall::mult_hop_omp(Field& w5, const void* f5) const{
  int nid = omp_get_num_threads();
  int tid = omp_get_thread_num();

  BGQThread_Barrier(0,nid);
  if(tid == 0) BGWilson_DW_Init(prms_.N5_,prms_.mq_,prms_.M0_, 
				(double*)&prms_.dp_[0],(double*)&prms_.dm_[0],
				(double*)&prms_.bs_[0],(double*)&prms_.cs_[0],
				(double*)&prms_.es_[0],(double*)&prms_.fs_[0]);
  BGQThread_Barrier(0,nid);
  double kappa = 0.5/(4.0+prms_.M0_);
  double* u = const_cast<Field *>(Dw_->getGaugeField_ptr())->getaddr(0);
  BGWilson_DW_Mult_hop(w5.getaddr(0),(void*)u,(void*)f5,kappa,BGWILSON_DIRAC);
}

void Dirac_DomainWall::mult_hop_dag_omp(Field& w5, const void* f5) const{
  int nid = omp_get_num_threads();
  int tid = omp_get_thread_num();

  BGQThread_Barrier(0,nid);

  if(tid == 0) BGWilson_DW_Init(prms_.N5_,prms_.mq_,prms_.M0_,
				(double*)&prms_.dp_[0],(double*)&prms_.dm_[0],
				(double*)&prms_.bs_[0],(double*)&prms_.cs_[0],
				(double*)&prms_.es_[0],(double*)&prms_.fs_[0]);
  BGQThread_Barrier(0,nid);

  double kappa = 0.5/(4.0+prms_.M0_);
  double* u = const_cast<Field *>(Dw_->getGaugeField_ptr())->getaddr(0);

  BGWilson_DW_Mult_hop_dag(w5.getaddr(0),(void*)u,(void*)f5,
			   kappa,BGWILSON_DIRAC);
}
