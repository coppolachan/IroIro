/*!--------------------------------------------------------------------------
 * @file dirac_DomainWall_BGQ.cpp
 *
 * @brief Definition of some BGQ optimized methods for Dirac_optimalDomainWall (5d operator)
 *
 *-------------------------------------------------------------------------*/
#include "dirac_DomainWall.hpp"
#include "Solver/solver.hpp"

#ifdef IBM_BGQ_WILSON
#include "bgqwilson.h"
#include "bgqthread.h"
#include <omp.h>
#include <hwi/include/bqc/A2_inlines.h>
#endif

typedef struct FermionSpinor{
  double _Complex v[12];
}Spinor;
typedef struct GaugeConf{
  double _Complex v[9];
}GaugePtr;


//#define ENABLE_THREADING

void Dirac_optimalDomainWall::mult_hop(Field& w5, const Field& f5) const{
  typedef struct FermionSpinor{
    double _Complex v[12];
  }Spinor;
  typedef struct GaugeConf{
    double _Complex v[9];
  }GaugePtr;

  Field lpf(f4size_);
  Field lmf(f4size_), ey(f4size_), fy(f4size_);
  Field temp(fsize_);

  Spinor* w_ptr    = (Spinor*)malloc(sizeof(double)*f4size_);
  Spinor* v_ptr    = (Spinor*)malloc(sizeof(double)*f4size_);

  Spinor* lpf_ptr  = (Spinor*)lpf.getaddr(0);
  Spinor* lmf_ptr  = (Spinor*)lmf.getaddr(0);
  Spinor* ey_ptr   = (Spinor*)ey.getaddr(0);
  Spinor* fy_ptr   = (Spinor*)fy.getaddr(0);
  Spinor* temp_ptr = (Spinor*)temp.getaddr(0);
  Spinor* temp_ptr_bdry = (Spinor*)temp.getaddr((N5_-1)*f4size_);
  Spinor* w5_ptr_bdry   = (Spinor*)w5.getaddr((N5_-1)*f4size_);

  GaugePtr* pU       = (GaugePtr*)(const_cast<Field *>(u_)->getaddr(0));
  Spinor* f5_ptr     = (Spinor*)(const_cast<Field&>(f5).getaddr(0));

  double mass_fact   = 4.0+M0_;

  double minus_kappa = -Dw_.getKappa();
  double doe_factor = minus_kappa*mass_fact;
  double fact= 1.0/(Params.dp_[N5_-1] +mq_*Params.dm_[N5_-2]*Params.es_[N5_-2]);

  register int Nvol = CommonPrms::instance()->Nvol()/2;

  Format::Format_F Fformat = Dw_.get_fermionFormat();

  int tid, nid;
  int is, ns;


  ///////////////////////////////  doe
#pragma omp parallel private(tid,nid,is,ns) firstprivate(f5_ptr, temp_ptr)
  {
    nid = omp_get_num_threads();
    tid = omp_get_thread_num();
    is = (Nvol * tid / nid);
    ns = Nvol/nid;
    
    // s = 0
    BGWilsonLA_Proj_P(lpf_ptr + is,const_cast<Field&>(f5).getaddr(Fformat.index(0,0,N5_-1))+is,ns);
    BGWilsonLA_Proj_M(lmf_ptr + is,const_cast<Field&>(f5).getaddr(Fformat.index(0,0,1))+is,ns);
    BGWilsonLA_MultScalar_Add(lpf_ptr + is,lmf_ptr + is,-mq_,ns);
    BGWilsonLA_AXPBY(v_ptr + is, f5_ptr + is, lpf_ptr + is, Params.bs_[0], Params.cs_[0], ns);

    //BGWilson_MultEO(w_ptr, pU, v_ptr, minus_kappa , 2, BGWILSON_DIRAC);
    //BGWilsonLA_MultScalar(temp_ptr+is, w_ptr+is, mass_fact, ns);
    BGWilson_MultEO(temp_ptr, pU, v_ptr, doe_factor , 2, BGWILSON_DIRAC);

  
    for(int s=1; s<N5_-1; ++s) {
      f5_ptr   += Nvol;
      temp_ptr += Nvol;
      
      BGWilsonLA_Proj_P(lpf_ptr+is,const_cast<Field&>(f5).getaddr(Fformat.index(0,0,s-1))+is,ns);
      BGWilsonLA_Proj_M(lmf_ptr+is,const_cast<Field&>(f5).getaddr(Fformat.index(0,0,s+1))+is,ns);
      BGWilsonLA_Add(lpf_ptr + is,lmf_ptr + is ,ns);
      BGWilsonLA_AXPBY(v_ptr + is , f5_ptr + is, lpf_ptr + is , Params.bs_[s],Params.cs_[s],ns);
     
      //BGWilson_MultEO(w_ptr, pU, v_ptr, minus_kappa , 2, BGWILSON_DIRAC);
      //BGWilsonLA_MultScalar(temp_ptr+is, w_ptr+is, mass_fact, ns);
      BGWilson_MultEO(temp_ptr, pU, v_ptr, doe_factor , 2, BGWILSON_DIRAC);
    }
  
    // s = N5-1
    f5_ptr   += Nvol;
    temp_ptr += Nvol;
    
    BGWilsonLA_Proj_P(lpf_ptr+is,const_cast<Field&>(f5).getaddr(Fformat.index(0,0,N5_-2))+is,ns);
    BGWilsonLA_Proj_M(lmf_ptr+is,const_cast<Field&>(f5).getaddr(Fformat.index(0,0,0))+is,ns);
    BGWilsonLA_MultAddScalar(lpf_ptr+is,lmf_ptr+is,-mq_,ns);  
    BGWilsonLA_AXPBY(v_ptr+is, f5_ptr+is, lpf_ptr+is, Params.bs_[N5_-1],Params.cs_[N5_-1],ns);

    //BGWilson_MultEO(w_ptr, pU, v_ptr, minus_kappa , 2, BGWILSON_DIRAC);
    //BGWilsonLA_MultScalar(temp_ptr+is, w_ptr+is, mass_fact, ns);
    //uint64_t time1 = GetTimeBase();
    //for(int i = 0; i < 100; ++i)
    BGWilson_MultEO(temp_ptr, pU, v_ptr, doe_factor , 2, BGWILSON_DIRAC); 
    //uint64_t time2 = GetTimeBase();
    //if (tid == 0){
    //  CCIO::cout << "Dw clocks : " << time2-time1 << "\n";
    //}
  


    /// 5d hopping term 

    //time1 = GetTimeBase();

    BGWilsonLA_Proj_M(ey_ptr+is,(Spinor*)temp.getaddr(0)+is,ns);
    BGWilsonLA_MultAddScalar(temp_ptr_bdry+is,     ey_ptr+is,-mq_* Params.es_[0],ns);

    for (int s=1; s<N5_-1; ++s) {
      Spinor* temp_ptr   = (Spinor*)temp.getaddr(s*f4size_);
      double fact_lpf = (Params.dm_[s]/Params.dp_[s-1]);
      double fact_ey =  mq_*Params.es_[s];
      
      BGWilsonLA_Proj_P(lpf_ptr+is ,(Spinor*)temp.getaddr(Fformat.index(0,0,s-1))+is ,ns);
      BGWilsonLA_Proj_M(ey_ptr+is  ,(Spinor*)temp.getaddr(Fformat.index(0,0,s))+is   ,ns);
      
      BGWilsonLA_MultAddScalar(temp_ptr+is     ,lpf_ptr+is, fact_lpf,ns);
      BGWilsonLA_MultAddScalar(temp_ptr_bdry+is,ey_ptr+is, -fact_ey,ns);

    }
    BGWilsonLA_Proj_P(lpf_ptr+is,(Spinor*)temp.getaddr(Fformat.index(0,0,N5_-2))+is,ns);
    BGWilsonLA_MultAddScalar(temp_ptr_bdry+is,     lpf_ptr+is,(Params.dm_[N5_-1]/Params.dp_[N5_-2]),ns);
    
    BGWilsonLA_MultScalar(temp_ptr_bdry+is, temp_ptr_bdry+is, fact, ns);
    
    for(int s=N5_-2; s>=0; --s) {
      Spinor* temp_ptr   = (Spinor*)temp.getaddr(s*f4size_);
      BGWilsonLA_Proj_M(lmf_ptr+is,(Spinor*)temp.getaddr(Fformat.index(0,0,s+1))+is,ns);
      BGWilsonLA_Proj_P(fy_ptr+is,(Spinor*)temp.getaddr(Fformat.index(0,0,N5_-1))+is,ns);
      BGWilsonLA_MultAddScalar(temp_ptr+is,     lmf_ptr+is,Params.dm_[s],ns);
      BGWilsonLA_MultAddScalar(temp_ptr+is,     fy_ptr+is,-mq_*Params.fs_[s],ns);
      BGWilsonLA_MultScalar(temp_ptr+is, temp_ptr+is, 1.0/ Params.dp_[s], ns);
    }

    /*
      time2 = GetTimeBase();
      if (tid == 0){
      CCIO::cout << "Dee clocks : " << time2-time1 << "\n";
      }
    */

    ////////////////deo
    for(int s=0; s<N5_; ++s) {
      Spinor* temp_ptr   = (Spinor*)temp.getaddr(s*f4size_);
      Spinor* w5_ptr   = (Spinor*)w5.getaddr(s*f4size_);
      BGWilsonLA_Proj_P(lpf_ptr+is,(Spinor*)temp.getaddr(Fformat.index(0,0,(s+N5_-1)%N5_))+is,ns);
      BGWilsonLA_Proj_M(lmf_ptr+is,(Spinor*)temp.getaddr(Fformat.index(0,0,(s+1)%N5_))+is,ns);
      
      if(s==0){
	BGWilsonLA_MultScalar_Add(lpf_ptr+is,lmf_ptr+is,-mq_,ns);
      }
      else if(s==N5_-1){
	BGWilsonLA_MultAddScalar(lpf_ptr+is,lmf_ptr+is,-mq_,ns);
      }
      else{
	BGWilsonLA_Add(lpf_ptr+is,lmf_ptr+is,ns);
      }
      BGWilsonLA_AXPBY(v_ptr+is, temp_ptr+is, lpf_ptr+is,
		       Params.bs_[s],Params.cs_[s],ns);
      
      BGWilson_MultEO(w_ptr, pU, v_ptr, minus_kappa , 1, BGWILSON_DIRAC);
      
      BGWilsonLA_MultScalar(w5_ptr+is, w_ptr+is, mass_fact, ns);
    }

    BGWilsonLA_Proj_M(ey_ptr+is,(Spinor*)w5.getaddr(Fformat.index(0,0,0))+is,ns);
    BGWilsonLA_MultAddScalar(w5_ptr_bdry+is,     ey_ptr+is,-mq_* Params.es_[0],ns);
    
    for (int s=1; s<N5_-1; ++s) {
      Spinor* w5_ptr   = (Spinor*)w5.getaddr(s*f4size_);
      double fact_lpf = (Params.dm_[s]/Params.dp_[s-1]);
      double fact_ey =  mq_*Params.es_[s];
      
      BGWilsonLA_Proj_P(lpf_ptr+is,(Spinor*)w5.getaddr(Fformat.index(0,0,s-1))+is,ns);
      BGWilsonLA_Proj_M(ey_ptr+is,(Spinor*)w5.getaddr(Fformat.index(0,0,s))+is,ns);
      
      BGWilsonLA_MultAddScalar(w5_ptr+is,     lpf_ptr+is,fact_lpf,ns);
      BGWilsonLA_MultAddScalar(w5_ptr_bdry+is,ey_ptr+is,-fact_ey,ns);
      
    }


      
    BGWilsonLA_Proj_P(lpf_ptr+is,(Spinor*)w5.getaddr(Fformat.index(0,0,N5_-2))+is,ns);
    BGWilsonLA_MultAddScalar(w5_ptr_bdry+is,     lpf_ptr+is,(Params.dm_[N5_-1]/Params.dp_[N5_-2]),ns);
    BGWilsonLA_MultScalar(w5_ptr_bdry+is, w5_ptr_bdry+is, fact, ns);
    
    for(int s=N5_-2; s>=0; --s) {
      Spinor* w5_ptr   = (Spinor*)w5.getaddr(s*f4size_);
      BGWilsonLA_Proj_M(lmf_ptr+is  ,(Spinor*)w5.getaddr(Fformat.index(0,0,s+1))+is,ns);
      BGWilsonLA_Proj_P(fy_ptr+is   ,(Spinor*)w5.getaddr(Fformat.index(0,0,N5_-1))+is,ns);
      BGWilsonLA_MultAddScalar(w5_ptr+is,     lmf_ptr+is,Params.dm_[s],ns);
      BGWilsonLA_MultAddScalar(w5_ptr+is,     fy_ptr+is,-mq_*Params.fs_[s],ns);
      BGWilsonLA_MultScalar(w5_ptr+is, w5_ptr+is, 1.0/ Params.dp_[s], ns);
    }
     
    //BGWilsonLA_MultScalar_Add((Spinor*)(w5.getaddr(0))+is*N5_,(Spinor*)(const_cast<Field&>(f5).getaddr(0))+is*N5_, -1.0, ns*N5_);
  
    Spinor* f5_ptr   = (Spinor*)(const_cast<Field&>(f5).getaddr(0));
    Spinor* w5_ptr   = (Spinor*)w5.getaddr(0);
    for(int s=0; s<N5_; ++s) {
      BGWilsonLA_MultScalar_Add(w5_ptr+is,f5_ptr+is, -1.0, ns);
      w5_ptr += Nvol;
      f5_ptr += Nvol;
    }    
  
  }  

  free(v_ptr);
  free(w_ptr);

}
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_THREADING
void Dirac_optimalDomainWall::mult_hop_omp(Field& w5, const void* f5) const{

  //Field temp(fsize_);
  register int Nvol = CommonPrms::instance()->Nvol()/2;
  int tid, nid;
  int is, ns;
  nid = omp_get_num_threads();
  tid = omp_get_thread_num();
  is = (Nvol * tid / nid);
  ns = Nvol/nid;

  int sp_size = f4size_;
  /*  
  Spinor* w_ptr    = (Spinor*)malloc(sizeof(double)*sp_size);
  Spinor* v_ptr    = (Spinor*)malloc(sizeof(double)*sp_size);
  Spinor* lpf_ptr  = (Spinor*)malloc(sizeof(double)*sp_size);
  Spinor* lmf_ptr  = (Spinor*)malloc(sizeof(double)*sp_size);
  Spinor* ey_ptr   = (Spinor*)malloc(sizeof(double)*sp_size);
  Spinor* fy_ptr   = (Spinor*)malloc(sizeof(double)*sp_size);
  */
  
  Spinor* w_ptr    = (Spinor*)BGQThread_Malloc(sizeof(double)*sp_size, nid);
  Spinor* v_ptr    = (Spinor*)BGQThread_Malloc(sizeof(double)*sp_size, nid);
  Spinor* lpf_ptr  = (Spinor*)BGQThread_Malloc(sizeof(double)*sp_size, nid);
  Spinor* lmf_ptr  = (Spinor*)BGQThread_Malloc(sizeof(double)*sp_size, nid);
  Spinor* ey_ptr   = (Spinor*)BGQThread_Malloc(sizeof(double)*sp_size, nid);
  Spinor* fy_ptr   = (Spinor*)BGQThread_Malloc(sizeof(double)*sp_size, nid);
  
  //Spinor* temp_ptr_base = (Spinor*)temp.getaddr(0);
  Spinor* temp_ptr_base = (Spinor*)BGQThread_Malloc(sizeof(double)*fsize_,nid);

  BGQThread_Barrier(0,nid);

  Spinor* temp_ptr = temp_ptr_base;
  Spinor* temp_ptr_bdry = temp_ptr_base+(N5_-1)*(Nvol);
  Spinor* w5_ptr_bdry   = (Spinor*)w5.getaddr((N5_-1)*sp_size);

  GaugePtr* pU       = (GaugePtr*)(const_cast<Field *>(u_)->getaddr(0));
  Spinor* f5_ptr     = (Spinor*)f5;

  double mass_fact   = 4.0+M0_;

  double minus_kappa = -Dw_.getKappa();
  double doe_factor = minus_kappa*mass_fact;
  double fact= 1.0/(Params.dp_[N5_-1] +mq_*Params.dm_[N5_-2]*Params.es_[N5_-2]);


  Format::Format_F Fformat = Dw_.get_fermionFormat();

  
  ///////////////////////////////  doe
  // s = 0
  BGWilsonLA_Proj_P(lpf_ptr + is,(f5_ptr+Fformat.index(0,0,N5_-1))+is,ns);
  BGWilsonLA_Proj_M(lmf_ptr + is,(f5_ptr+Fformat.index(0,0,1))+is    ,ns);
  BGWilsonLA_MultScalar_Add(lpf_ptr + is,lmf_ptr + is,-mq_,ns);
  BGWilsonLA_AXPBY(v_ptr + is, f5_ptr + is, lpf_ptr + is, Params.bs_[0], Params.cs_[0], ns);
  
  BGWilson_MultEO(temp_ptr, pU, v_ptr, doe_factor , 2, BGWILSON_DIRAC);
    
  for(int s=1; s<N5_-1; ++s) {
    f5_ptr   += Nvol;
    temp_ptr += Nvol;
      
    BGWilsonLA_Proj_P(lpf_ptr+is,(f5_ptr+Fformat.index(0,0,s-1))+is,ns);
    BGWilsonLA_Proj_M(lmf_ptr+is,(f5_ptr+Fformat.index(0,0,s+1))+is,ns);
    BGWilsonLA_Add(lpf_ptr + is,lmf_ptr + is ,ns);
    BGWilsonLA_AXPBY(v_ptr + is , f5_ptr + is, lpf_ptr + is , Params.bs_[s],Params.cs_[s],ns);
     
    BGWilson_MultEO(temp_ptr, pU, v_ptr, doe_factor , 2, BGWILSON_DIRAC);
  }
  
  // s = N5-1
  f5_ptr   += Nvol;
  temp_ptr += Nvol;
    
  BGWilsonLA_Proj_P(lpf_ptr+is,(f5_ptr+Fformat.index(0,0,N5_-2))+is,ns);
  BGWilsonLA_Proj_M(lmf_ptr+is,(f5_ptr+Fformat.index(0,0,0))+is    ,ns);
  BGWilsonLA_MultAddScalar(lpf_ptr+is,lmf_ptr+is,-mq_,ns);  
  BGWilsonLA_AXPBY(v_ptr+is, f5_ptr+is, lpf_ptr+is, Params.bs_[N5_-1],Params.cs_[N5_-1],ns);

  BGWilson_MultEO(temp_ptr, pU, v_ptr, doe_factor , 2, BGWILSON_DIRAC); 
  ///////////////////////////
  /// 5d hopping term 
  BGWilsonLA_Proj_M(ey_ptr+is,temp_ptr_base+is,ns);
  BGWilsonLA_MultAddScalar(temp_ptr_bdry+is,     ey_ptr+is,-mq_* Params.es_[0],ns);

  for (int s=1; s<N5_-1; ++s) {
    Spinor* temp_ptr   = temp_ptr_base+s*Nvol;//(Spinor*)temp.getaddr(s*f4size_);
    double fact_lpf = (Params.dm_[s]/Params.dp_[s-1]);
    double fact_ey =  mq_*Params.es_[s];
      
    BGWilsonLA_Proj_P(lpf_ptr+is ,temp_ptr_base+(s-1)*Nvol+is ,ns);
    BGWilsonLA_Proj_M(ey_ptr+is  ,temp_ptr_base+s*Nvol+is   ,ns);
      
    BGWilsonLA_MultAddScalar(temp_ptr+is     ,lpf_ptr+is, fact_lpf,ns);
    BGWilsonLA_MultAddScalar(temp_ptr_bdry+is,ey_ptr+is, -fact_ey,ns);

  }
  BGWilsonLA_Proj_P(lpf_ptr+is,temp_ptr_base+(N5_-2)*Nvol+is,ns);
  BGWilsonLA_MultAddScalar(temp_ptr_bdry+is,     lpf_ptr+is,(Params.dm_[N5_-1]/Params.dp_[N5_-2]),ns);
    
  BGWilsonLA_MultScalar(temp_ptr_bdry+is, temp_ptr_bdry+is, fact, ns);
    
  for(int s=N5_-2; s>=0; --s) {
    Spinor* temp_ptr   = temp_ptr_base+s*Nvol;
    BGWilsonLA_Proj_M(lmf_ptr+is,temp_ptr_base+(s+1)*Nvol+is,ns);
    BGWilsonLA_Proj_P(fy_ptr+is,temp_ptr_base+(N5_-1)*Nvol+is,ns);
    BGWilsonLA_MultAddScalar(temp_ptr+is,     lmf_ptr+is,Params.dm_[s],ns);
    BGWilsonLA_MultAddScalar(temp_ptr+is,     fy_ptr+is,-mq_*Params.fs_[s],ns);
    BGWilsonLA_MultScalar(temp_ptr+is, temp_ptr+is, 1.0/ Params.dp_[s], ns);
  }

  ////////////////deo
  for(int s=0; s<N5_; ++s) {
    Spinor* temp_ptr   = temp_ptr_base+s*Nvol;
    Spinor* w5_ptr   = (Spinor*)w5.getaddr(s*f4size_);
    BGWilsonLA_Proj_P(lpf_ptr+is,temp_ptr_base+((s+N5_-1)%N5_)*Nvol+is,ns);
    BGWilsonLA_Proj_M(lmf_ptr+is,temp_ptr_base+((s+1)%N5_)*Nvol+is,ns);
      
    if(s==0){
      BGWilsonLA_MultScalar_Add(lpf_ptr+is,lmf_ptr+is,-mq_,ns);
    }
    else if(s==N5_-1){
      BGWilsonLA_MultAddScalar(lpf_ptr+is,lmf_ptr+is,-mq_,ns);
    }
    else{
      BGWilsonLA_Add(lpf_ptr+is,lmf_ptr+is,ns);
    }
    BGWilsonLA_AXPBY(v_ptr+is, temp_ptr+is, lpf_ptr+is,
		     Params.bs_[s],Params.cs_[s],ns);
      
    BGWilson_MultEO(w_ptr, pU, v_ptr, minus_kappa , 1, BGWILSON_DIRAC);
      
    BGWilsonLA_MultScalar(w5_ptr+is, w_ptr+is, mass_fact, ns);
  }

  BGWilsonLA_Proj_M(ey_ptr+is,(Spinor*)w5.getaddr(Fformat.index(0,0,0))+is,ns);
  BGWilsonLA_MultAddScalar(w5_ptr_bdry+is,     ey_ptr+is,-mq_* Params.es_[0],ns);
    
  for (int s=1; s<N5_-1; ++s) {
    Spinor* w5_ptr   = (Spinor*)w5.getaddr(s*f4size_);
    double fact_lpf = (Params.dm_[s]/Params.dp_[s-1]);
    double fact_ey =  mq_*Params.es_[s];
      
    BGWilsonLA_Proj_P(lpf_ptr+is,(Spinor*)w5.getaddr(Fformat.index(0,0,s-1))+is,ns);
    BGWilsonLA_Proj_M(ey_ptr+is,(Spinor*)w5.getaddr(Fformat.index(0,0,s))+is,ns);
      
    BGWilsonLA_MultAddScalar(w5_ptr+is,     lpf_ptr+is,fact_lpf,ns);
    BGWilsonLA_MultAddScalar(w5_ptr_bdry+is,ey_ptr+is,-fact_ey,ns);
      
  }
      
  BGWilsonLA_Proj_P(lpf_ptr+is,(Spinor*)w5.getaddr(Fformat.index(0,0,N5_-2))+is,ns);
  BGWilsonLA_MultAddScalar(w5_ptr_bdry+is,     lpf_ptr+is,(Params.dm_[N5_-1]/Params.dp_[N5_-2]),ns);
  BGWilsonLA_MultScalar(w5_ptr_bdry+is, w5_ptr_bdry+is, fact, ns);
    
  for(int s=N5_-2; s>=0; --s) {
    Spinor* w5_ptr   = (Spinor*)w5.getaddr(s*f4size_);
    BGWilsonLA_Proj_M(lmf_ptr+is  ,(Spinor*)w5.getaddr(Fformat.index(0,0,s+1))+is,ns);
    BGWilsonLA_Proj_P(fy_ptr+is   ,(Spinor*)w5.getaddr(Fformat.index(0,0,N5_-1))+is,ns);
    BGWilsonLA_MultAddScalar(w5_ptr+is,     lmf_ptr+is,Params.dm_[s],ns);
    BGWilsonLA_MultAddScalar(w5_ptr+is,     fy_ptr+is,-mq_*Params.fs_[s],ns);
    BGWilsonLA_MultScalar(w5_ptr+is, w5_ptr+is, 1.0/ Params.dp_[s], ns);
  }
     
  //BGWilsonLA_MultScalar_Add((Spinor*)(w5.getaddr(0))+is*N5_,(Spinor*)(const_cast<Field&>(f5).getaddr(0))+is*N5_, -1.0, ns*N5_);
  
  f5_ptr   = (Spinor*)f5;
  Spinor* w5_ptr   = (Spinor*)w5.getaddr(0);
  for(int s=0; s<N5_; ++s) {
    BGWilsonLA_MultScalar_Add(w5_ptr+is,f5_ptr+is, -1.0, ns);
    w5_ptr += Nvol;
    f5_ptr += Nvol;
  }    

  BGQThread_Barrier(0, nid);

#pragma omp single
  {
    free(v_ptr);
    free(w_ptr);
    free(lpf_ptr);
    free(lmf_ptr);
    free(ey_ptr);
    free(fy_ptr);
    free(temp_ptr_base);
  }


}
#endif


void Dirac_optimalDomainWall::mult_hop_dag(Field& w5, const Field& f5) const{
  typedef struct FermionSpinor{
    double _Complex v[12];
  }Spinor;
  typedef struct GaugeConf{
    double _Complex v[9];
  }GaugePtr;

  register int Nvol = CommonPrms::instance()->Nvol()/2;
  Field lpf(f4size_), ey(f4size_), lmf(f4size_), v(f4size_), w(f4size_);
  Field temp(f5),v5(fsize_);

  Spinor* w5_ptr  = (Spinor*)w5.getaddr(0);
  Spinor* temp_ptr = (Spinor*)temp.getaddr(0);
  Spinor* v_ptr   = (Spinor*)v.getaddr(0);
  Spinor* lpf_ptr = (Spinor*)lpf.getaddr(0);
  Spinor* lmf_ptr = (Spinor*)lmf.getaddr(0);
  Spinor* ey_ptr  = (Spinor*)ey.getaddr(0);
  Spinor* w_ptr = (Spinor*)w.getaddr(0);

  int spin_idx;
  double cs, bs;
  double minus_kappa = -Dw_.getKappa();
  double* pU       = const_cast<Field *>(u_)->getaddr(0);


  int tid, nid;
  int is, ns;


  ///////////////////////////////  doe
#pragma omp parallel private(tid,nid,is,ns) firstprivate(temp_ptr)
  {
    nid = omp_get_num_threads();
    tid = omp_get_thread_num();
    is = (Nvol * tid / nid);
    ns = Nvol/nid;
    
    // 5d hopping term
    BGWilsonLA_MultScalar(temp_ptr+is, temp_ptr+is, 1.0/ Params.dp_[0], ns);
    
    for(int s=1; s<N5_-1; ++s){
      Spinor* temp_ptr   = (Spinor*)temp.getaddr(s*f4size_);
      BGWilsonLA_Proj_M(lmf_ptr+is,
			(Spinor*)(const_cast<Field&>(temp).getaddr(Dw_.get_fermionFormat().index(0,0,s-1)))+is,
			ns);
      BGWilsonLA_MultAddScalar(temp_ptr+is,     lmf_ptr+is,Params.dm_[s-1],ns);
      BGWilsonLA_MultScalar(temp_ptr+is, temp_ptr+is, 1.0/ Params.dp_[s], ns);   
    } 
    
    temp_ptr   += (N5_-1)*Nvol;
    BGWilsonLA_Proj_M(v_ptr+is,
		      (Spinor*)(const_cast<Field&>(temp).getaddr(Dw_.get_fermionFormat().index(0,0,N5_-2)))+is,ns);   
    BGWilsonLA_MultAddScalar(temp_ptr+is, v_ptr+is,Params.dm_[N5_-2],ns);
    for(int s=0; s<N5_-1; ++s) {
      BGWilsonLA_Proj_P(ey_ptr+is,
			(Spinor*)(const_cast<Field&>(temp).getaddr(Dw_.get_fermionFormat().index(0,0,s)))+is,ns); 
      BGWilsonLA_MultAddScalar(temp_ptr+is, ey_ptr+is,-mq_*Params.fs_[s],ns);
    }
    BGWilsonLA_MultScalar(temp_ptr+is, temp_ptr+is,1.0/(Params.dp_[N5_-1] +mq_*Params.dm_[N5_-2]*Params.es_[N5_-2]),ns);
    
    
    for(int s=N5_-2; s>=0; --s){ 
      Spinor* temp_ptr   = (Spinor*)temp.getaddr(s*f4size_);
      double fact_lpf = (Params.dm_[s+1]/Params.dp_[s]);
      
      BGWilsonLA_Proj_P(lpf_ptr+is,
			(Spinor*)(const_cast<Field&>(temp).getaddr(Dw_.get_fermionFormat().index(0,0,s+1)))+is,ns); 
      BGWilsonLA_Proj_M(ey_ptr+is,
			(Spinor*)(const_cast<Field&>(temp).getaddr(Dw_.get_fermionFormat().index(0,0,N5_-1)))+is,ns);
      BGWilsonLA_AXPBYPZ(temp_ptr+is, lpf_ptr+is,ey_ptr+is, temp_ptr+is,
			 fact_lpf, - mq_*Params.es_[s],ns);
    }

    //mult_off_diag deo
    for(int s=0; s<N5_; ++s){
      spin_idx = s*f4size_;
      Spinor* temp_ptr = (Spinor*)temp.getaddr(spin_idx);
      Spinor* w5_ptr = (Spinor*)w5.getaddr(spin_idx);
      Spinor* v5_ptr = (Spinor*)v5.getaddr(spin_idx);
      
      bs = (4.0+M0_)*Params.bs_[s];
      cs = (4.0+M0_)*Params.cs_[s];
      
      BGWilson_MultEO_Dag(w_ptr, pU, temp_ptr, minus_kappa , 2, BGWILSON_DIRAC);
      
      BGWilsonLA_MultScalar(w5_ptr+is, w_ptr+is, bs,ns);
      BGWilsonLA_MultScalar(v5_ptr+is, w_ptr+is, cs,ns);
    }


    for(int s = 0; s < N5_; ++s){
      spin_idx = s*f4size_;
      Spinor* w5_ptr = (Spinor*)w5.getaddr(spin_idx);
      Spinor* v5_ptr = (Spinor*)v5.getaddr(spin_idx);

      BGWilsonLA_Proj_P(lpf_ptr+is,
			(Spinor*)(const_cast<Field&>(v5).getaddr(Dw_.get_fermionFormat().index(0,0,(s+1)%N5_)))+is,
			ns);
      BGWilsonLA_Proj_M(lmf_ptr+is,
			(Spinor*)(const_cast<Field&>(v5).getaddr(Dw_.get_fermionFormat().index(0,0,(s+N5_-1)%N5_)))+is,
			ns);
    
      if(s == N5_-1){
	BGWilsonLA_AXPBYPZ(w5_ptr+is, lpf_ptr+is,lmf_ptr+is,
			   w5_ptr+is, -mq_,1.0,ns);
      }
      else if(s==0){
	BGWilsonLA_AXPBYPZ(w5_ptr+is, lpf_ptr+is,lmf_ptr+is, w5_ptr+is,
			   1.0,-mq_,ns);
      }
      else{
	BGWilsonLA_AXPBYPZ(w5_ptr+is, lpf_ptr+is,lmf_ptr+is, w5_ptr+is,
			   1.0,1.0,ns);
      }
    }

    // 5d hopping term
    BGWilsonLA_MultScalar(w5_ptr+is, w5_ptr+is, 1.0/ Params.dp_[0], ns);

    for(int s=1; s<N5_-1; ++s){
      w5_ptr   = (Spinor*)w5.getaddr(s*f4size_);
      BGWilsonLA_Proj_M(lmf_ptr+is,
			(Spinor*)(const_cast<Field&>(w5).getaddr(Dw_.get_fermionFormat().index(0,0,s-1)))+is,ns);
      BGWilsonLA_MultAddScalar(w5_ptr+is,     lmf_ptr+is,Params.dm_[s-1],ns);
      BGWilsonLA_MultScalar(w5_ptr+is, w5_ptr+is, 1.0/ Params.dp_[s], ns);   
    } 

    w5_ptr   = (Spinor*)w5.getaddr((N5_-1)*f4size_); 
    BGWilsonLA_Proj_M(v_ptr+is,
		      (Spinor*)(const_cast<Field&>(w5).getaddr(Dw_.get_fermionFormat().index(0,0,N5_-2)))+is,ns);   
    BGWilsonLA_MultAddScalar(w5_ptr+is, v_ptr+is,Params.dm_[N5_-2],ns);
    for(int s=0; s<N5_-1; ++s) {
      BGWilsonLA_Proj_P(ey_ptr+is,
			(Spinor*)(const_cast<Field&>(w5).getaddr(Dw_.get_fermionFormat().index(0,0,s)))+is,ns); 
      BGWilsonLA_MultAddScalar(w5_ptr+is, ey_ptr+is,-mq_*Params.fs_[s],ns);
    }
    BGWilsonLA_MultScalar(w5_ptr+is, w5_ptr+is,1.0/(Params.dp_[N5_-1] +mq_*Params.dm_[N5_-2]*Params.es_[N5_-2]),ns);


    for(int s=N5_-2; s>=0; --s){ 
      Spinor* w5_ptr   = (Spinor*)w5.getaddr(s*f4size_);
      double fact_lpf = (Params.dm_[s+1]/Params.dp_[s]);
    
      BGWilsonLA_Proj_P(lpf_ptr+is,
			(Spinor*)(const_cast<Field&>(w5).getaddr(Dw_.get_fermionFormat().index(0,0,s+1)))+is,ns); 
      BGWilsonLA_Proj_M(ey_ptr+is,
			(Spinor*)(const_cast<Field&>(w5).getaddr(Dw_.get_fermionFormat().index(0,0,N5_-1)))+is,ns);
      BGWilsonLA_AXPBYPZ(w5_ptr+is, lpf_ptr+is,ey_ptr+is, w5_ptr+is,
			 fact_lpf, - mq_*Params.es_[s],ns);
    }
  
    //mult_off_diag doe
    for(int s=0; s<N5_; ++s){
      spin_idx = s*f4size_;
      Spinor* w5_ptr = (Spinor*)w5.getaddr(spin_idx);
      Spinor* v5_ptr = (Spinor*)v5.getaddr(spin_idx);
    
      bs = (4.0+M0_)*Params.bs_[s];
      cs = (4.0+M0_)*Params.cs_[s];
    
      BGWilson_MultEO_Dag(w_ptr, pU, w5_ptr, minus_kappa , 1, BGWILSON_DIRAC);
    
      BGWilsonLA_MultScalar(w5_ptr+is, w_ptr+is, bs,ns);
      BGWilsonLA_MultScalar(v5_ptr+is, w_ptr+is, cs,ns);
    
    }

    for(int s = 0; s < N5_; ++s){
      spin_idx = s*f4size_;
      Spinor* w5_ptr = (Spinor*)w5.getaddr(spin_idx);
      Spinor* v5_ptr = (Spinor*)v5.getaddr(spin_idx);

      BGWilsonLA_Proj_P(lpf_ptr+is,
			(Spinor*)(const_cast<Field&>(v5).getaddr(Dw_.get_fermionFormat().index(0,0,(s+1)%N5_)))+is,
			ns);
      BGWilsonLA_Proj_M(lmf_ptr+is,
			(Spinor*)(const_cast<Field&>(v5).getaddr(Dw_.get_fermionFormat().index(0,0,(s+N5_-1)%N5_)))+is,
			ns);
    
      if(s == N5_-1){
	BGWilsonLA_AXPBYPZ(w5_ptr+is, lpf_ptr+is,lmf_ptr+is,
			   w5_ptr+is, -mq_,1.0,ns);
      }
      else if(s==0){
	BGWilsonLA_AXPBYPZ(w5_ptr+is, lpf_ptr+is,lmf_ptr+is, w5_ptr+is,
			   1.0,-mq_,ns);
      }
      else{
	BGWilsonLA_AXPBYPZ(w5_ptr+is, lpf_ptr+is,lmf_ptr+is, w5_ptr+is,
			   1.0,1.0,ns);
      }
    }

    Spinor* f5_ptr   = (Spinor*)(const_cast<Field&>(f5).getaddr(0));
    Spinor* w5_ptr   = (Spinor*)w5.getaddr(0);
    for(int s=0; s<N5_; ++s) {
      BGWilsonLA_MultScalar_Add(w5_ptr+is,f5_ptr+is, -1.0, ns);
      w5_ptr += Nvol;
      f5_ptr += Nvol;
    }    
  
  }

  

}
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_THREADING
void Dirac_optimalDomainWall::mult_hop_dag_omp(Field& w5, const Field& f5) const{
  register int Nvol = CommonPrms::instance()->Nvol()/2;

  Field temp(f5);
  Field v5(fsize_);

  // 4d vectors
  Spinor* w_ptr    = (Spinor*)malloc(sizeof(double)*f4size_);
  Spinor* v_ptr    = (Spinor*)malloc(sizeof(double)*f4size_);
  Spinor* lpf_ptr  = (Spinor*)malloc(sizeof(double)*f4size_);
  Spinor* lmf_ptr  = (Spinor*)malloc(sizeof(double)*f4size_);
  Spinor* ey_ptr   = (Spinor*)malloc(sizeof(double)*f4size_);

  // 5d vector
  Spinor* v5_ptr_base =  (Spinor*)v5.getaddr(0);
  //Spinor* v5_ptr_base = (Spinor*)malloc(sizeof(double)*fsize_);

  Spinor* w5_ptr  = (Spinor*)w5.getaddr(0);
  Spinor* temp_ptr_base = (Spinor*)temp.getaddr(0);
  Spinor* temp_ptr = temp_ptr_base;
  Spinor* f5_ptr   = (Spinor*)(const_cast<Field&>(f5).getaddr(0));
  
  int spin_idx;
  double cs, bs;
  double minus_kappa = -Dw_.getKappa();
  double* pU       = const_cast<Field *>(u_)->getaddr(0);

  int tid, nid;
  int is, ns;


  ///////////////////////////////  doe
  nid = omp_get_num_threads();
  tid = omp_get_thread_num();
  is = (Nvol * tid / nid);
  ns = Nvol/nid;
    
  // 5d hopping term
  BGWilsonLA_MultScalar(temp_ptr+is, temp_ptr+is, 1.0/ Params.dp_[0], ns);
    
  for(int s=1; s<N5_-1; ++s){
    Spinor* temp_ptr   = temp_ptr_base+s*Nvol;
    BGWilsonLA_Proj_M(lmf_ptr+is,
		      temp_ptr_base+(s-1)*Nvol+is,
		      ns);
    BGWilsonLA_MultAddScalar(temp_ptr+is,     lmf_ptr+is,Params.dm_[s-1],ns);
    BGWilsonLA_MultScalar(temp_ptr+is, temp_ptr+is, 1.0/ Params.dp_[s], ns);   
  } 
    
  temp_ptr   += (N5_-1)*Nvol;
  BGWilsonLA_Proj_M(v_ptr+is,
		    temp_ptr_base+(N5_-2)*Nvol+is,
		    ns);   
  BGWilsonLA_MultAddScalar(temp_ptr+is, v_ptr+is,Params.dm_[N5_-2],ns);
  for(int s=0; s<N5_-1; ++s) {
    BGWilsonLA_Proj_P(ey_ptr+is,
		      temp_ptr_base+(s)*Nvol+is,
		      ns); 
    BGWilsonLA_MultAddScalar(temp_ptr+is, ey_ptr+is,-mq_*Params.fs_[s],ns);
  }
  BGWilsonLA_MultScalar(temp_ptr+is, temp_ptr+is,1.0/(Params.dp_[N5_-1] +mq_*Params.dm_[N5_-2]*Params.es_[N5_-2]),ns);
    
    
  for(int s=N5_-2; s>=0; --s){ 
    Spinor* temp_ptr   = temp_ptr_base+s*Nvol;
    double fact_lpf = (Params.dm_[s+1]/Params.dp_[s]);
      
    BGWilsonLA_Proj_P(lpf_ptr+is,
		      temp_ptr_base+(s+1)*Nvol+is,
		      ns); 
    BGWilsonLA_Proj_M(ey_ptr+is,
		      temp_ptr_base+(N5_-1)*Nvol+is,
		      ns);
    BGWilsonLA_AXPBYPZ(temp_ptr+is, lpf_ptr+is,ey_ptr+is, temp_ptr+is,
		       fact_lpf, - mq_*Params.es_[s],ns);
  }

  //mult_off_diag deo
  for(int s=0; s<N5_; ++s){
    spin_idx = s*f4size_;
    Spinor* temp_ptr = temp_ptr_base+s*Nvol;
    Spinor* w5_ptr = (Spinor*)w5.getaddr(spin_idx);
    Spinor* v5_ptr = v5_ptr_base+s*Nvol;//(Spinor*)v5.getaddr(spin_idx);
      
    bs = (4.0+M0_)*Params.bs_[s];
    cs = (4.0+M0_)*Params.cs_[s];
      
    BGWilson_MultEO_Dag(w_ptr, pU, temp_ptr, minus_kappa , 2, BGWILSON_DIRAC);
      
    BGWilsonLA_MultScalar(w5_ptr+is, w_ptr+is, bs,ns);
    BGWilsonLA_MultScalar(v5_ptr+is, w_ptr+is, cs,ns);
  }


  for(int s = 0; s < N5_; ++s){
    spin_idx = s*f4size_;
    Spinor* w5_ptr = (Spinor*)w5.getaddr(spin_idx);
    Spinor* v5_ptr = v5_ptr_base+s*Nvol;//(Spinor*)v5.getaddr(spin_idx);

    BGWilsonLA_Proj_P(lpf_ptr+is,
		      v5_ptr_base+((s+1)%N5_)*Nvol+is,
		      //(Spinor*)(const_cast<Field&>(v5).getaddr(Dw_.get_fermionFormat().index(0,0,(s+1)%N5_)))+is,
		      ns);
    BGWilsonLA_Proj_M(lmf_ptr+is,
		      v5_ptr_base+((s+N5_-1)%N5_)*Nvol+is,
		      //(Spinor*)(const_cast<Field&>(v5).getaddr(Dw_.get_fermionFormat().index(0,0,(s+N5_-1)%N5_)))+is,
		      ns);
    
    if(s == N5_-1){
      BGWilsonLA_AXPBYPZ(w5_ptr+is, lpf_ptr+is,lmf_ptr+is,
			 w5_ptr+is, -mq_,1.0,ns);
    }
    else if(s==0){
      BGWilsonLA_AXPBYPZ(w5_ptr+is, lpf_ptr+is,lmf_ptr+is, w5_ptr+is,
			 1.0,-mq_,ns);
    }
    else{
      BGWilsonLA_AXPBYPZ(w5_ptr+is, lpf_ptr+is,lmf_ptr+is, w5_ptr+is,
			 1.0,1.0,ns);
    }
  }

  // 5d hopping term
  BGWilsonLA_MultScalar(w5_ptr+is, w5_ptr+is, 1.0/ Params.dp_[0], ns);

  for(int s=1; s<N5_-1; ++s){
    w5_ptr   = (Spinor*)w5.getaddr(s*f4size_);
    BGWilsonLA_Proj_M(lmf_ptr+is,
		      (Spinor*)(const_cast<Field&>(w5).getaddr(Dw_.get_fermionFormat().index(0,0,s-1)))+is,ns);
    BGWilsonLA_MultAddScalar(w5_ptr+is,     lmf_ptr+is,Params.dm_[s-1],ns);
    BGWilsonLA_MultScalar(w5_ptr+is, w5_ptr+is, 1.0/ Params.dp_[s], ns);   
  } 

  w5_ptr   = (Spinor*)w5.getaddr((N5_-1)*f4size_); 
  BGWilsonLA_Proj_M(v_ptr+is,
		    (Spinor*)(const_cast<Field&>(w5).getaddr(Dw_.get_fermionFormat().index(0,0,N5_-2)))+is,ns);   
  BGWilsonLA_MultAddScalar(w5_ptr+is, v_ptr+is,Params.dm_[N5_-2],ns);
  for(int s=0; s<N5_-1; ++s) {
    BGWilsonLA_Proj_P(ey_ptr+is,
		      (Spinor*)(const_cast<Field&>(w5).getaddr(Dw_.get_fermionFormat().index(0,0,s)))+is,ns); 
    BGWilsonLA_MultAddScalar(w5_ptr+is, ey_ptr+is,-mq_*Params.fs_[s],ns);
  }
  BGWilsonLA_MultScalar(w5_ptr+is, w5_ptr+is,1.0/(Params.dp_[N5_-1] +mq_*Params.dm_[N5_-2]*Params.es_[N5_-2]),ns);


  for(int s=N5_-2; s>=0; --s){ 
    Spinor* w5_ptr   = (Spinor*)w5.getaddr(s*f4size_);
    double fact_lpf = (Params.dm_[s+1]/Params.dp_[s]);
    
    BGWilsonLA_Proj_P(lpf_ptr+is,
		      (Spinor*)(const_cast<Field&>(w5).getaddr(Dw_.get_fermionFormat().index(0,0,s+1)))+is,ns); 
    BGWilsonLA_Proj_M(ey_ptr+is,
		      (Spinor*)(const_cast<Field&>(w5).getaddr(Dw_.get_fermionFormat().index(0,0,N5_-1)))+is,ns);
    BGWilsonLA_AXPBYPZ(w5_ptr+is, lpf_ptr+is,ey_ptr+is, w5_ptr+is,
		       fact_lpf, - mq_*Params.es_[s],ns);
  }
  
  //mult_off_diag doe
  for(int s=0; s<N5_; ++s){
    spin_idx = s*f4size_;
    Spinor* w5_ptr = (Spinor*)w5.getaddr(spin_idx);
    Spinor* v5_ptr = (Spinor*)v5.getaddr(spin_idx);
    
    bs = (4.0+M0_)*Params.bs_[s];
    cs = (4.0+M0_)*Params.cs_[s];
    
    BGWilson_MultEO_Dag(w_ptr, pU, w5_ptr, minus_kappa , 1, BGWILSON_DIRAC);
    
    BGWilsonLA_MultScalar(w5_ptr+is, w_ptr+is, bs,ns);
    BGWilsonLA_MultScalar(v5_ptr+is, w_ptr+is, cs,ns);
    
  }

  for(int s = 0; s < N5_; ++s){
    spin_idx = s*f4size_;
    Spinor* w5_ptr = (Spinor*)w5.getaddr(spin_idx);
    Spinor* v5_ptr = (Spinor*)v5.getaddr(spin_idx);

    BGWilsonLA_Proj_P(lpf_ptr+is,
		      (Spinor*)(const_cast<Field&>(v5).getaddr(Dw_.get_fermionFormat().index(0,0,(s+1)%N5_)))+is,
		      ns);
    BGWilsonLA_Proj_M(lmf_ptr+is,
		      (Spinor*)(const_cast<Field&>(v5).getaddr(Dw_.get_fermionFormat().index(0,0,(s+N5_-1)%N5_)))+is,
		      ns);
    
    if(s == N5_-1){
      BGWilsonLA_AXPBYPZ(w5_ptr+is, lpf_ptr+is,lmf_ptr+is,
			 w5_ptr+is, -mq_,1.0,ns);
    }
    else if(s==0){
      BGWilsonLA_AXPBYPZ(w5_ptr+is, lpf_ptr+is,lmf_ptr+is, w5_ptr+is,
			 1.0,-mq_,ns);
    }
    else{
      BGWilsonLA_AXPBYPZ(w5_ptr+is, lpf_ptr+is,lmf_ptr+is, w5_ptr+is,
			 1.0,1.0,ns);
    }
  }


  w5_ptr   = (Spinor*)w5.getaddr(0);
  for(int s=0; s<N5_; ++s) {
    BGWilsonLA_MultScalar_Add(w5_ptr+is,f5_ptr+is, -1.0, ns);
    w5_ptr += Nvol;
    f5_ptr += Nvol;
  }    
  
  
  free(w_ptr);
  free(v_ptr);
  free(lpf_ptr);
  free(lmf_ptr);
  free(ey_ptr);
  

}
#endif



#ifdef IBM_BGQ_WILSON
// CG solver optimized for BGQ
void Dirac_optimalDomainWall::solve_eo_5d(Field& w5, const Field& b, SolverOutput& Out, int MaxIter, double GoalPrecision) const{
  
  //#if VERBOSITY>=SOLV_ITER_VERB_LEVEL
  CCIO::header("CG_BGQ solver start");
  //#endif

#ifdef ENABLE_THREADING


  BGQThread_Init(); //initializing BGQ fast threading routines


  Out.Msg = "CG solver";
  Out.Iterations = -1;
  double kernel_timing = 0.0; 
  int tid, nid;
  int is, ns;
  register int Nvol = CommonPrms::instance()->Nvol()/2;
  Field temp(fsize_), temp2(fsize_), s(fsize_); //eventually eliminated
  TIMING_START;
  
  Field x = b;//initial condition
  Field r = b;//initial residual


  double* x_ptr = x.getaddr(0);
  double* r_ptr = r.getaddr(0);
  double* s_ptr = s.getaddr(0);
  double* temp_ptr = temp.getaddr(0);
  #pragma omp parallel 
  { 
  mult_hop_omp(temp,x_ptr);
  }
  mult_hop_dag_omp(temp2,temp);
  //}

  r -= temp2;

  Field p = r;
  double rr = r*r;// (r,r)
  double snorm = b.norm();
  snorm = 1.0/snorm;

  double* p_ptr = p.getaddr(0);
  int v_size = x.size()/24; // << assumes 24 elements

  //#if VERBOSITY>1
  CCIO::cout<<" Snorm = "<< snorm << std::endl;
  CCIO::cout<<" Init  = "<< rr*snorm<< std::endl;
  //#endif

  nid = 1;//omp_get_num_threads();
  tid = 0;//omp_get_thread_num();
  is = (v_size * tid / nid);
  ns = v_size/nid;


  for(int it = 0; it < MaxIter; ++it){
    mult_hop_omp(temp,p_ptr);
    mult_hop_dag_omp(s,temp);

    ///////////////////////////////////////////////////
 
    double pap;
    BGWilsonLA_DotProd(&pap,p_ptr,s_ptr,v_size);
    pap = Communicator::instance()->reduce_sum(pap);
    double rrp = rr;
    double cr = rrp/pap;// (r,r)/(p,Ap)

    //  x += cr*p; // x = x + cr * p
    BGWilsonLA_MultAddScalar(x_ptr,p_ptr,cr,v_size);
    //  r -= cr*s; // r_k = r_k - cr * Ap
    BGWilsonLA_MultAddScalar(r_ptr,s_ptr,-cr,v_size);
  
    //  rr = r*r; // rr = (r_k,r_k)
    BGWilsonLA_Norm(&rr,r_ptr,v_size);
    rr = Communicator::instance()->reduce_sum(rr);
    //  p *= rr/rrp; // p = p*(r_k,r_k)/(r,r)
    //  p += r; // p = p + p*(r_k,r_k)/(r,r)
    BGWilsonLA_MultScalar_Add(p_ptr,r_ptr,rr/rrp,v_size);
    /////////////////////////////////////////////////////
    //#if VERBOSITY>1
    CCIO::cout<< std::setw(5)<< "["<<it<<"] "
	      << std::setw(20) << rr*snorm<< "\n";
    //#endif    
    if(rr*snorm < GoalPrecision){
      Out.Iterations = it;
      break;
    }
  }

  

  if(Out.Iterations == -1) {
    CCIO::cout<<" Not converged. Current residual: "<< rr*snorm << "\n";
    abort();
  }

  //p = opr_->mult(x);
  mult_hop(temp,x);
  mult_hop_dag(p,temp);

  p -= b;
  Out.diff = p.norm();

  w5 = x;

  TIMING_END(Out.timing);
  CCIO::cout << "Kernel section timing: "<< kernel_timing << "\n";

#endif

}
#endif
