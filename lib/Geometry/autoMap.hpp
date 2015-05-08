/*!
  @file autoMap.hpp
  @brief Declares AutoMap class for the shift
*/
#ifndef MAPPING_H_
#define MAPPING_H_

#include "include/common_fields.hpp"
#include "include/errors.hpp"
#include "Communicator/communicator.hpp"
#include "siteMap.hpp"
#include "myMap.hpp"

#include <vector>
#include <valarray>

#ifdef IBM_BGQ_WILSON
#include <omp.h>
#include "bgqthread.h"

#include <hwi/include/bqc/A2_inlines.h>
#endif


//   ____________   ->-    ____________
//  |*|          |   mu   |          |*|
//  |*|  bulk_b_ |        | bulk_t_  |*|
//  |*|__________|        |__________|*|
//   ^ bdry_b_                bdry_t_ ^

namespace Mapping{
  static double *class_send_f, *class_recv_f;
  static double *class_send_b, *class_recv_b;


  struct SpatialTag{};

  struct EOtag{};
  struct OEtag{};

  struct Forward{};
  struct Backward{};

  class AutoMap{
    int is_initialized_;
    int dir_;
    int Nvol_;
    std::vector<int> bdry_t_;
    std::vector<int> bulk_t_;
    std::vector<int> bdry_b_;
    std::vector<int> bulk_b_;

    AutoMap(){}//default constructor

  public:
    // using the regular site indexing <-
    AutoMap(int dir):bdry_t_(SiteMap::shiftSite.bdry_map(dir,Top)),
		     bulk_t_(SiteMap::shiftSite.bulk_map(dir,Top)),
		     bdry_b_(SiteMap::shiftSite.bdry_map(dir,Btm)),
		     bulk_b_(SiteMap::shiftSite.bulk_map(dir,Btm)),
		     dir_(dir),
		     Nvol_(CommonPrms::instance()->Nvol()),
		     is_initialized_(true){
      if(bdry_t_.empty()||bulk_t_.empty()
	 ||bdry_b_.empty()||bulk_b_.empty()) is_initialized_=false;
    }
    template<typename MAP>
    AutoMap(int dir,const MAP& mm):bdry_t_(mm.bdry_map(dir,Top)),
				   bulk_t_(mm.bulk_map(dir,Top)),
				   bdry_b_(mm.bdry_map(dir,Btm)),
				   bulk_b_(mm.bulk_map(dir,Btm)),
				   dir_(dir),
				   is_initialized_(true){
      if(bdry_t_.empty()||bulk_t_.empty()
	 ||bdry_b_.empty()||bulk_b_.empty()) is_initialized_=false;
    }
    /*
    template<typename FIELD>
    FIELD operator()(const FIELD& Fin,Forward)const{
      if (!is_initialized_){
	ErrorString msg;
	msg << "The maps is not initialized correctly. dir=" <<dir_<<"\n";
	Errors::BaseErr("Initialization missing", msg);
      }
      FIELD Fout(Fin.Nvol());
      std::valarray<double> recv_bdry(bdry_t_.size()*Fin.Nin()*Fin.Nex());
      Communicator::instance()->transfer_fw(recv_bdry,
					    Fin.data[Fin.get_sub(bdry_b_)],
					    dir_);
      
      Fout.data.set(Fin.get_sub(bdry_t_),recv_bdry);
      Fout.data.set(Fin.get_sub(bulk_t_),Fin.data[Fin.get_sub(bulk_b_)]);
      return Fout;
    }
    */
    template<typename FIELD>
    FIELD operator()(const FIELD& Fin,Forward)const{
      FIELD Fout(Fin.Nvol());
      operator()(Fout,Fin,Forward());
      return Fout;
    }

    template<typename FIELD>
    void operator()(FIELD& Fout,const FIELD& Fin,Forward)const{
      if(!is_initialized_){
	ErrorString msg;
	msg << "The maps is not initialized correctly. dir=" <<dir_<<"\n";
	Errors::BaseErr("Initialization missing", msg);
      }
      std::valarray<double> recv_bdry(bdry_t_.size()*Fin.Nin()*Fin.Nex());
      Communicator::instance()->transfer_fw(recv_bdry,
					    Fin.data[Fin.get_sub(bdry_b_)],
					    dir_);
      Fout.data.set(Fin.get_sub(bdry_t_),recv_bdry);
      Fout.data.set(Fin.get_sub(bulk_t_),Fin.data[Fin.get_sub(bulk_b_)]);
    }
    /*
    template<typename FIELD>
    FIELD operator()(const FIELD& Fin,Backward) const{
      if (!is_initialized_){
	ErrorString msg;
	msg << "The maps is not initialized correctly. dir=" <<dir_<<"\n";
	Errors::BaseErr("Initialization missing", msg);
      }
      FIELD Fout(Fin.Nvol());   
      std::valarray<double> recv_bdry(bdry_b_.size()*Fin.Nin()*Fin.Nex());
      
      Communicator::instance()->transfer_bk(recv_bdry,
					    Fin.data[Fin.get_sub(bdry_t_)],
					    dir_);

      Fout.data.set(Fin.get_sub(bdry_b_),recv_bdry);
      Fout.data.set(Fin.get_sub(bulk_b_),Fin.data[Fin.get_sub(bulk_t_)]);
      return Fout;
    }
    */

    template<typename FIELD>
    FIELD operator()(const FIELD& Fin,Backward)const{
      FIELD Fout(Fin.Nvol());
      operator()(Fout,Fin,Backward());
      return Fout;
    }

    template<typename FIELD>
    void operator()(FIELD& Fout,const FIELD& Fin,Backward) const{
      if (!is_initialized_){
	ErrorString msg;
	msg << "The maps is not initialized correctly. dir=" <<dir_<<"\n";
	Errors::BaseErr("Initialization missing", msg);
      }
      std::valarray<double> recv_bdry(bdry_b_.size()*Fin.Nin()*Fin.Nex());
      Communicator::instance()->transfer_bk(recv_bdry,
					    Fin.data[Fin.get_sub(bdry_t_)],
					    dir_);
      Fout.data.set(Fin.get_sub(bdry_b_),recv_bdry);
      Fout.data.set(Fin.get_sub(bulk_b_),Fin.data[Fin.get_sub(bulk_t_)]);
    }

    /// for experimental use
    void operator()(GaugeField1D& Fout,const GaugeField1D& Fin,Forward)const;
    void operator()(GaugeField1D& Fout,const GaugeField1D& Fin,Backward)const;
    void operator()(GaugeField1D& Fout,const GaugeField& Fin,int mu,Forward)const;
    void operator()(GaugeField1D& Fout,const GaugeField& Fin,int mu,Backward)const;

    // for improvements assuming the memory alignment
    // using OpenMP and BGQ threading
    template <typename CommonField1D>
    void operator()(CommonField1D& Fout,const double* Fin,Forward)const{
#ifdef IBM_BGQ_WILSON
    if(!omp_in_parallel()){
#endif
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
#ifdef IBM_BGQ_WILSON
    }else{
      //variables declared here are private by default           
      int tID = omp_get_thread_num();
      int nID = omp_get_num_threads();
      int Nin = Fout.Nin();
      int block   = bdry_t_.size()/nID;
      int bulk_bl = bulk_b_.size()/nID;

#pragma omp master
      {      
	class_send_f = 
	  (double*)malloc(bdry_t_.size()*Nin*sizeof(double));
	class_recv_f =
	  (double*)malloc(bdry_t_.size()*Nin*sizeof(double));
      }
      //BGQThread_Barrier(0, nID);      
#pragma omp barrier 
      
      for(int b=0; b<block; ++b)
	for(int i=0; i<Nin; ++i)
	  class_send_f[(b+tID*block)*Nin+i] = Fin[bdry_b_[b+tID*block]*Nin+i];
      
      //BGQThread_Barrier(0, nID);
#pragma omp barrier 

      if(tID == 0) {
	Communicator::instance()->transfer_fw(class_recv_f,class_send_f,
					      bdry_t_.size()*Nin,dir_);
	Communicator::instance()->sync();
      }
      //BGQThread_Barrier(0, nID);
#pragma omp barrier 
      
      for(int b=0; b<block; ++b)
	for(int i=0; i<Nin; ++i)
	  Fout.data.set(Fout.format.index(i,bdry_t_[b+tID*block]),
			class_recv_f[(b+tID*block)*Nin+i]);
      	
      for(int b=0; b<bulk_bl; ++b)
	for(int i=0; i<Nin; ++i)
	  Fout.data.set(Fout.format.index(i,bulk_t_[b+tID*bulk_bl]),
			Fin[bulk_b_[b+tID*bulk_bl]*Nin+i]);
      	
      //BGQThread_Barrier(0, nID);
#pragma omp barrier 
#pragma omp master
      {   
    	free(class_send_f);
	free(class_recv_f);

      }
#pragma omp barrier
      

    }
#endif


    };


    // Only BGQ
   template <typename CommonField1D>
   void operator()(CommonField1D& Fout,const double* Fin,Backward)const{
#ifdef IBM_BGQ_WILSON
    if(!omp_in_parallel()){
#endif
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
#ifdef IBM_BGQ_WILSON      
    }else{
      //variables declared here are private by default                           
      int tID = omp_get_thread_num();
      int nID = omp_get_num_threads();
      int Nin = Fout.Nin();
      int block   = bdry_b_.size()/nID;
      int bulk_bl = bulk_t_.size()/nID;
      
      //double *class_send_b, *class_recv_b;
  
#pragma omp master
      {  

     class_send_b =
	(double*)malloc(bdry_b_.size()*Nin*sizeof(double));
      class_recv_b =
	(double*)malloc(bdry_b_.size()*Nin*sizeof(double));
      }
      //BGQThread_Barrier(0, nID);      
#pragma omp barrier    

      for(int b=0; b<block; ++b)
	for(int i=0; i<Nin; ++i)
	  class_send_b[(b+tID*block)*Nin+i] = Fin[bdry_t_[b+tID*block]*Nin+i];
      
      //BGQThread_Barrier(0, nID);
#pragma omp barrier    

      if(tID == 0) {
	Communicator::instance()->transfer_bk(class_recv_b,class_send_b,
					      bdry_b_.size()*Nin,dir_);
	Communicator::instance()->sync();
      }
      //BGQThread_Barrier(0, nID);
#pragma omp barrier    
      
      for(int b=0; b<block; ++b)
	for(int i=0; i<Nin; ++i)
	  Fout.data.set(Fout.format.index(i,bdry_b_[b+tID*block]),
			class_recv_b[(b+tID*block)*Nin+i]);
      
      for(int b=0; b<bulk_bl; ++b)
	for(int i=0; i<Nin; ++i)
	  Fout.data.set(Fout.format.index(i,bulk_b_[b+tID*bulk_bl]),
			Fin[bulk_t_[b+tID*bulk_bl]*Nin+i]);
      

      //BGQThread_Barrier(0, nID);
#pragma omp barrier    
#pragma omp master
      {

	free(class_send_b);
	free(class_recv_b);

      }
      //BGQThread_Barrier(0, nID);
#pragma omp barrier    
      
    }
#endif

   };
 

  };

  class AutoMap_EvenOdd{
    int is_initialized_;
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
       dir_(dir),
       is_initialized_(true){}

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
       dir_(dir),
       is_initialized_(true){}

    template<typename FIELD>
    FIELD operator()(const FIELD& Fin,Forward)const{
      if (!is_initialized_){
	ErrorString msg;
	msg << "The maps have not been initialized yet\n";
	Errors::BaseErr("Initialization missing", msg);
      }
      FIELD Fout(Fin.Nvol());
      std::valarray<double> recv_bdry(recv_bdry_t_.size()*Fin.Nin()*Fin.Nex());

      Communicator::instance()->transfer_fw(recv_bdry,
					    Fin.data[Fin.get_sub(send_bdry_b_)],
					    dir_);
      Fout.data.set(Fin.get_sub(recv_bdry_t_),recv_bdry);
      Fout.data.set(Fin.get_sub(recv_bulk_t_),Fin.data[Fin.get_sub(send_bulk_b_)]);
      return Fout;
    }

    template<typename FIELD>
    FIELD operator()(const FIELD& Fin,Backward)const{
      if (!is_initialized_){
	ErrorString msg;
	msg << "The maps have not been initialized yet\n";
	Errors::BaseErr("Initialization missing", msg);
      }
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

    /// for improvements assuming the memory alinement
    void operator()(GaugeField1D& Fout,const double* Fin,Forward)const;
    void operator()(GaugeField1D& Fout,const double* Fin,Backward)const;
    // now a hack for FermionField generalize to all fields
    void operator()(FermionField& Fout,const double* Fin,Forward)const;



  };
}
#endif

