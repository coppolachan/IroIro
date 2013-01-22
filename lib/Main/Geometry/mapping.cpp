/*!
  @file mapping.cpp
  @brief Definition of ShiftField class
*/
#include "shiftField.hpp"
#include "include/macros.hpp"

#ifdef IBM_BGQ_WILSON
#include <omp.h>
#include "bgqthread.h"
#endif


namespace Mapping{

#ifdef IBM_BGQ_WILSON
  //////// for BGQ ///////
  void AutoMap::operator()(GaugeField1D& Fout,const double* Fin,Forward)const{
    if (!omp_in_parallel()){
      int Nin = Fout.Nin();
      int bdsize = Nin*bdry_t_.size();
      double send_bdry[bdsize],recv_bdry[bdsize];
      
      for(int b=0; b<bdry_t_.size(); ++b)
	for(int i=0; i<Nin; ++i)
	  send_bdry[b*Nin+i] = Fin[bdry_b_[b]*Nin+i];
      
      Communicator::instance()->transfer_fw(recv_bdry,send_bdry,bdsize,dir_);
      Communicator::instance()->sync();
	
	for(int b=0; b<bdry_t_.size(); ++b)      
	  for(int i=0; i<Nin; ++i)
	    Fout.data.set(Fout.format.index(i,bdry_t_[b]),recv_bdry[b*Nin+i]);
	
	for(int b=0; b<bulk_b_.size(); ++b)      
	  for(int i=0; i<Nin; ++i)
	    Fout.data.set(Fout.format.index(i,bulk_t_[b]),Fin[bulk_b_[b]*Nin+i]);
      }
      else
	{
	  //variables declared here are private by default
	  int tID, nID;
	  tID = omp_get_thread_num();
	  nID = omp_get_num_threads();
	  int Nin = Fout.Nin();
	  int block   = bdry_t_.size()/nID;
	  int bulk_bl = bulk_b_.size()/nID;
	  int bdsize = Nin*bdry_t_.size()/nID;
	  double* class_send = 
	    (double*)BGQThread_Malloc(bdry_t_.size()*Nin*sizeof(double), nID);
	  double* class_recv =
	    (double*)BGQThread_Malloc(bdry_t_.size()*Nin*sizeof(double), nID);

	  for(int b=0; b<block; ++b)
	    for(int i=0; i<Nin; ++i){
	      class_send[(b+tID*block)*Nin+i] = Fin[bdry_b_[b+tID*block]*Nin+i];
	    }
	  
	  BGQThread_Barrier(0, nID);

	  if(tID == 0)
	    Communicator::instance()->transfer_fw(class_recv,class_send,bdry_t_.size()*Nin,dir_);
	  
	  BGQThread_Barrier(0, nID);

	  for(int b=0; b<block; ++b)      
	    for(int i=0; i<Nin; ++i)
	      Fout.data.set(Fout.format.index(i,bdry_t_[b+tID*block]),class_recv[(b+tID*block)*Nin+i]);

	  for(int b=0; b<bulk_bl; ++b)      
	    for(int i=0; i<Nin; ++i)
	      Fout.data.set(Fout.format.index(i,bulk_t_[b+tID*bulk_bl]),Fin[bulk_b_[b+tID*bulk_bl]*Nin+i]);
	  
	  BGQThread_Barrier(0, nID);
	  BGQThread_Free(class_send, tID);
	  BGQThread_Free(class_recv, tID);


	}
    }
  #endif

#ifdef IBM_BGQ_WILSON
  // assumes internal indexing for the Fields
  void AutoMap::operator()(GaugeField1D& Fout,const double* Fin,Backward)const{
    if (!omp_in_parallel()){
      int Nin = Fout.Nin();
      int bdsize = Nin*bdry_b_.size();
      double send_bdry[bdsize],recv_bdry[bdsize];
      
      for(int b=0; b<bdry_b_.size(); ++b)
        for(int i=0; i<Nin; ++i)
          send_bdry[b*Nin+i] = Fin[bdry_t_[b]*Nin+i];
      
      Communicator::instance()->transfer_bk(recv_bdry,send_bdry,bdsize,dir_);
      Communicator::instance()->sync();
      
      for(int b=0; b<bdry_b_.size(); ++b)      
        for(int i=0; i<Nin; ++i)
          Fout.data.set(Fout.format.index(i,bdry_b_[b]),recv_bdry[b*Nin+i]);
      
      for(int b=0; b<bulk_t_.size(); ++b)      
        for(int i=0; i<Nin; ++i)
          Fout.data.set(Fout.format.index(i,bulk_b_[b]),Fin[bulk_t_[b]*Nin+i]);
    }
    else 
      {
	//variables declared here are private by default
	int tID, nID;
	tID = omp_get_thread_num();
	nID = omp_get_num_threads();
	int Nin = Fout.Nin();
	int block   = bdry_b_.size()/nID;
	int bulk_bl = bulk_t_.size()/nID;
	double* class_send =
	  (double*)BGQThread_Malloc(bdry_b_.size()*Nin*sizeof(double), nID);
	double* class_recv = 
	  (double*)BGQThread_Malloc(bdry_b_.size()*Nin*sizeof(double), nID);
	 
	 for(int b=0; b<block; ++b)
	   for(int i=0; i<Nin; ++i){
	     class_send[(b+tID*block)*Nin+i] = Fin[bdry_t_[b+tID*block]*Nin+i];
	   }
	 
	 BGQThread_Barrier(0, nID);
	 
	 if(tID == 0)
	   Communicator::instance()->transfer_bk(class_recv,class_send,bdry_b_.size()*Nin,dir_);
	 
	 BGQThread_Barrier(0, nID);
	 
	 for(int b=0; b<block; ++b)      
	   for(int i=0; i<Nin; ++i)
	     Fout.data.set(Fout.format.index(i,bdry_b_[b+tID*block]),class_recv[(b+tID*block)*Nin+i]);
	 
	 for(int b=0; b<bulk_bl; ++b)      
	   for(int i=0; i<Nin; ++i)
	     Fout.data.set(Fout.format.index(i,bulk_b_[b+tID*bulk_bl]),Fin[bulk_t_[b+tID*bulk_bl]*Nin+i]);
	 
	 BGQThread_Barrier(0, nID);
	 BGQThread_Free(class_send, tID);
	 BGQThread_Free(class_recv, tID);
	 
      }
  }
#endif


  void ShiftField::init_maps(){
    if(maps_.size()==0){
      for(int dir=0; dir<NDIM_; ++dir)	
	maps_.push_back(AutoMap(dir));
    }
  }

  void ShiftField_eo::init_maps(){
    if(maps_.size()==0){
      for(int dir=0; dir<NDIM_; ++dir)	
	maps_.push_back(AutoMap_EvenOdd(dir,EOtag()));
    }
  }

  void ShiftField_oe::init_maps(){
    if(maps_.size()==0){
      for(int dir=0; dir<NDIM_; ++dir)	
	maps_.push_back(AutoMap_EvenOdd(dir,OEtag()));
    }
  }

  // global instance
  ShiftField shiftField;
  void init_shiftField(){ shiftField.init_maps();}

  ShiftField_eo shiftField_eo;
  ShiftField_oe shiftField_oe;

  void init_shiftField_EvenOdd(){ 
    shiftField_eo.init_maps();
    shiftField_oe.init_maps();
  }

}

