/*!
  @file action_gauge_wilson_force.bgq.cpp
  @brief Specialization of the md_force method for the ActionGaugeWilson class

  This is the BGQ optimized version with threading

  Time-stamp: <2013-12-13 11:10:41 cossu>
*/
#include "Action/action_gauge_wilson.hpp"
#include "Tools/fieldUtils.hpp"
#include "include/messages_macros.hpp"
#include <omp.h>
#include "bgqthread.h"

GaugeField ActionGaugeWilson::md_force(){
  using namespace Mapping;

  SUNmat pl;
  GaugeField force;
  GaugeField1D tmp; 
 
  GaugeField1D c;
  GaugeField1D WupMu, VupNu;

  double* VupNu_ptr = VupNu.data.getaddr(0);
  double* WupMu_ptr = WupMu.data.getaddr(0);
  double* c_ptr     = c.data.getaddr(0);
  double* tmp_ptr   = tmp.data.getaddr(0);
  double* U_ptr     = u_->data.getaddr(0);  
  double* U_mu_ptr; 
  double* U_nu_ptr; 
  double *inmat, *outmat;
  double trace;
  double factor = 0.5*Params.beta/NC_;

  BGQThread_Init();
  
#pragma omp parallel 
  {
    const int nid = omp_get_num_threads();
    const int tid = omp_get_thread_num();
    const int is = tid*Nvol_ / nid;
    const int ie = (tid + 1)*Nvol_ / nid;
    const int ns = ie - is;
    const int CC = NC_*NC_;
    const int str = is*CC;
    const int CC2 = 2*NC_*NC_;
    const int str2 = is*CC2;
    
    for(int mu = 0; mu < NDIM_; ++mu){
      BGWilsonSU3_MatZero(tmp_ptr+str2, ns);
      U_mu_ptr = U_ptr+mu*Nvol_*CC2; // assumes contiguity of U_mu
      for(int nu=0; nu< NDIM_; ++nu){
	if(nu == mu) continue;
	U_nu_ptr = U_ptr+nu*Nvol_*CC2; 
	
	shiftField(WupMu,U_nu_ptr,mu,Forward());
	shiftField(VupNu,U_mu_ptr,nu,Forward());
	
	BGWilsonSU3_MatMult_NND(c_ptr+str2, U_nu_ptr+str2, VupNu_ptr+str2, WupMu_ptr+str2, ns);
	BGWilsonSU3_MatMult_DNN(VupNu_ptr+str2, U_nu_ptr+str2, U_mu_ptr+str2, WupMu_ptr+str2, ns);
	shiftField(WupMu,VupNu_ptr,nu,Backward());
	BGWilsonSU3_MatAdd(tmp_ptr+str2,c_ptr+str2,ns);
	BGWilsonSU3_MatAdd(tmp_ptr+str2,WupMu_ptr+str2,ns);
      }
      BGWilsonSU3_MatMult_ND(c_ptr+str2,U_mu_ptr+str2,tmp_ptr+str2,ns);
      
      for(int site = is; site < (is+ns); ++site){
	inmat  = c_ptr + 2*NC_*NC_*site;//CC2*(site+mu*Nvol);
	outmat = force.data.getaddr(0)  +2*NC_*NC_*(site+mu*Nvol_);
	
	trace = factor*(inmat[1]+inmat[9]+inmat[17])/NC_;// imtrace
	for(int a=0; a<NC_; ++a){
	  for(int b=a; b<NC_; ++b){
	    
	    int ab = 2*(NC_*a+b);
	    int ba = 2*(NC_*b+a);
	    
	    *(outmat+ab)   = factor*0.5*(*(inmat+ab)  -*(inmat+ba));
	    *(outmat+ab+1) = factor*0.5*(*(inmat+ab+1)+*(inmat+ba+1));
	    
	    *(outmat+ba)   = -*(outmat+ab);
	    *(outmat+ba+1) = *(outmat+ab+1);
	    
	  }
	}
	outmat[1] -= trace;
	outmat[9] -= trace;
	outmat[17] -= trace;
	
      }     
    }
  }

  _MonitorMsg(ACTION_VERB_LEVEL, Action, force, "ActionGaugeWilson");

  return force;
}


