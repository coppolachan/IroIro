/*!
  @file mapping.hpp
  @brief Declares AutoMap class for the shift
*/
#ifndef MAPPING_H_
#define MAPPING_H_

#include "include/common_fields.hpp"
#include "Communicator/communicator.h"

#include "siteMap.hpp"

#include <vector>
#include <valarray>

#include <omp.h>
#include "bgqthread.h"

//   ____________   ->-    ____________
//  |*|          |   mu   |          |*|
//  |*|  bulk_b_ |        | bulk_t_  |*|
//  |*|__________|        |__________|*|
//   ^ bdry_b_                bdry_t_ ^

namespace Mapping{

  struct EOtag{};
  struct OEtag{};

  struct Forward{};
  struct Backward{};

  class AutoMap{
    int dir_;
    std::vector<int> bdry_t_;
    std::vector<int> bulk_t_;
    std::vector<int> bdry_b_;
    std::vector<int> bulk_b_;

    AutoMap(){}//default constructor

  public:
    // using the regular site indexing <-
    AutoMap(int dir)
      :bdry_t_(SiteMap::shiftSite.bdry_map(dir,Top)),
       bulk_t_(SiteMap::shiftSite.bulk_map(dir,Top)),
       bdry_b_(SiteMap::shiftSite.bdry_map(dir,Btm)),
       bulk_b_(SiteMap::shiftSite.bulk_map(dir,Btm)),
       dir_(dir){}

    template<typename FIELD>
    FIELD operator()(const FIELD& Fin,Forward)const{
      FIELD Fout(Fin.Nvol());
      std::valarray<double> recv_bdry(bdry_t_.size()*Fin.Nin()*Fin.Nex());
      
      Communicator::instance()->transfer_fw(recv_bdry,
					    Fin.data[Fin.get_sub(bdry_b_)],
					    dir_);
      Fout.data.set(Fin.get_sub(bdry_t_),recv_bdry);
      Fout.data.set(Fin.get_sub(bulk_t_),Fin.data[Fin.get_sub(bulk_b_)]);
      
      return Fout;
    }

    template<typename FIELD>
    FIELD operator()(const FIELD& Fin,Backward) const{
      FIELD Fout(Fin.Nvol());   
      std::valarray<double> recv_bdry(bdry_b_.size()*Fin.Nin()*Fin.Nex());
      
      Communicator::instance()->transfer_bk(recv_bdry,
					    Fin.data[Fin.get_sub(bdry_t_)],
					    dir_);

      Fout.data.set(Fin.get_sub(bdry_b_),recv_bdry);
      Fout.data.set(Fin.get_sub(bulk_b_),Fin.data[Fin.get_sub(bulk_t_)]);
     
      return Fout;
    }

    //////// for BGQ ///////
    void operator()(GaugeField1D& Fout,const double* Fin,Forward)const{
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
	  double* class_send = (double*)BGQThread_Malloc(bdry_t_.size()*Nin*sizeof(double), nID);
	  double* class_recv = (double*)BGQThread_Malloc(bdry_t_.size()*Nin*sizeof(double), nID);

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
  


    //for experimental use
    void operator()(GaugeField1D& Fout,const GaugeField1D& Fin,
		    int mu,Forward)const{
      std::valarray<double> recv_bdry(bdry_t_.size()*Fin.Nin());
      Communicator::instance()->transfer_fw(recv_bdry,
					    Fin.data[Fin.get_sub(bdry_b_)],
					    dir_);
      
      Fout.data.set(Fin.get_sub(bdry_t_),recv_bdry);
      Fout.data.set(Fin.get_sub(bulk_t_),Fin.data[Fin.get_sub(bulk_b_)]);
    }
    void operator()(GaugeField1D& Fout,const GaugeField& Fin,
		    int mu,Forward)const{
      std::valarray<double> recv_bdry(bdry_t_.size()*Fin.Nin());
      Communicator::instance()->transfer_fw(recv_bdry,
					    Fin.data[Fin.get_sub(bdry_b_,mu)],
					    dir_);
      
      Fout.data.set(Fin.get_sub(bdry_t_),recv_bdry);
      Fout.data.set(Fin.get_sub(bulk_t_),Fin.data[Fin.get_sub(bulk_b_,mu)]);
    }
    //for experimental use
    void operator()(GaugeField1D& Fout,const GaugeField1D& Fin,
		    int mu,Backward)const{
      std::valarray<double> recv_bdry(bdry_b_.size()*Fin.Nin());
      Communicator::instance()->transfer_bk(recv_bdry,
					    Fin.data[Fin.get_sub(bdry_t_)],
					    dir_);
      
      Fout.data.set(Fin.get_sub(bdry_b_),recv_bdry);
      Fout.data.set(Fin.get_sub(bulk_b_),Fin.data[Fin.get_sub(bulk_t_)]);
    }
    void operator()(GaugeField1D& Fout,const GaugeField& Fin,
		    int mu,Backward)const{
      std::valarray<double> recv_bdry(bdry_b_.size()*Fin.Nin());
      Communicator::instance()->transfer_bk(recv_bdry,
					    Fin.data[Fin.get_sub(bdry_t_,mu)],
					    dir_);
      
      Fout.data.set(Fin.get_sub(bdry_b_),recv_bdry);
      Fout.data.set(Fin.get_sub(bulk_b_),Fin.data[Fin.get_sub(bulk_t_,mu)]);
    }

    //////// for temporal use ///////
    // assumes internal indexing for the Fields
    void operator()(GaugeField1D& Fout,const double* Fin,Backward)const{
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
	 // variables declared here are private by default
	 int tID = omp_get_thread_num();
	 int nID = omp_get_num_threads();
	 int Nin = Fout.Nin();
	 int block   = bdry_b_.size()/nID;
	 int bulk_bl = bulk_t_.size()/nID;
	 // these will be shared
	 double* class_send = (double*)BGQThread_Malloc(bdry_b_.size()*Nin*sizeof(double), nID);
	 double* class_recv = (double*)BGQThread_Malloc(bdry_b_.size()*Nin*sizeof(double), nID);
	 
	 for(int b=0; b<block; ++b)
	   for(int i=0; i<Nin; ++i){
	     class_send[(b+tID*block)*Nin+i] = Fin[bdry_t_[b+tID*block]*Nin+i];
	   }
	 
	 BGQThread_Barrier(0, nID);
	 
	 // Transfer using the shared variables
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
  };

  class AutoMap_EvenOdd{
    int dir_;
    std::vector<int> recv_bdry_t_;
    std::vector<int> recv_bulk_t_;
    std::vector<int> recv_bdry_b_;
    std::vector<int> recv_bulk_b_;

    std::vector<int> send_bdry_t_;
    std::vector<int> send_bulk_t_;
    std::vector<int> send_bdry_b_;
    std::vector<int> send_bulk_b_;

    AutoMap_EvenOdd(){}//default constructor

  public:
    // using the e/o site indexing (e<-o)
    AutoMap_EvenOdd(int dir,EOtag)
      :send_bdry_t_(SiteMap::shiftSite_oe.bdry_map(dir,Top)),
       send_bulk_t_(SiteMap::shiftSite_oe.bulk_map(dir,Top)),
       send_bdry_b_(SiteMap::shiftSite_oe.bdry_map(dir,Btm)),
       send_bulk_b_(SiteMap::shiftSite_oe.bulk_map(dir,Btm)),
       recv_bdry_t_(SiteMap::shiftSite_eo.bdry_map(dir,Top)),
       recv_bulk_t_(SiteMap::shiftSite_eo.bulk_map(dir,Top)),
       recv_bdry_b_(SiteMap::shiftSite_eo.bdry_map(dir,Btm)),
       recv_bulk_b_(SiteMap::shiftSite_eo.bulk_map(dir,Btm)),
       dir_(dir){}

    // using the e/o site indexing (o<-e)
    AutoMap_EvenOdd(int dir,OEtag)
      :send_bdry_t_(SiteMap::shiftSite_eo.bdry_map(dir,Top)),
       send_bulk_t_(SiteMap::shiftSite_eo.bulk_map(dir,Top)),
       send_bdry_b_(SiteMap::shiftSite_eo.bdry_map(dir,Btm)),
       send_bulk_b_(SiteMap::shiftSite_eo.bulk_map(dir,Btm)),
       recv_bdry_t_(SiteMap::shiftSite_oe.bdry_map(dir,Top)),
       recv_bulk_t_(SiteMap::shiftSite_oe.bulk_map(dir,Top)),
       recv_bdry_b_(SiteMap::shiftSite_oe.bdry_map(dir,Btm)),
       recv_bulk_b_(SiteMap::shiftSite_oe.bulk_map(dir,Btm)),
       dir_(dir){}
    
    template<typename FIELD>
    FIELD operator()(const FIELD& Fin,Forward)const{
      FIELD Fout(Fin.Nvol());
      std::valarray<double> recv_bdry(recv_bdry_t_.size()*Fin.Nin()*Fin.Nex());

      Communicator::instance()->transfer_fw(recv_bdry,
					    Fin.data[Fin.get_sub(send_bdry_b_)],
					    dir_);
      Fout.data.set(Fin.get_sub(recv_bdry_t_),recv_bdry);
      Fout.data.set(Fin.get_sub(recv_bulk_t_),
		    Fin.data[Fin.get_sub(send_bulk_b_)]);
      return Fout;
    }

    template<typename FIELD>
    FIELD operator()(const FIELD& Fin,Backward)const{
      FIELD Fout(Fin.Nvol());
      std::valarray<double> recv_bdry(recv_bdry_b_.size()*Fin.Nin()*Fin.Nex());
      Communicator::instance()->transfer_bk(recv_bdry,
					    Fin.data[Fin.get_sub(send_bdry_t_)],
					    dir_);
      Fout.data.set(Fin.get_sub(recv_bdry_b_),recv_bdry);
      Fout.data.set(Fin.get_sub(recv_bulk_b_),
		    Fin.data[Fin.get_sub(send_bulk_t_)]);
      return Fout;
    }
  };
}
#endif

