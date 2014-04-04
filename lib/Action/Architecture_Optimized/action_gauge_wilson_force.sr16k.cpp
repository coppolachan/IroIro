/*!
  @file action_gauge_wilson_force.sr16k.cpp
  @brief Specialization of the md_force method for the ActionGaugeWilson class
  This is the SR16K optimized version
  Time-stamp: <2014-01-24 15:52:07 noaki>
*/
#include "Action/action_gauge_wilson.hpp"
#include "Tools/fieldUtils.hpp"
#include "include/messages_macros.hpp"
#include "srmwilson.h"

GaugeField ActionGaugeWilson::md_force(){
  using namespace FieldUtils;
  using namespace Mapping;

  GaugeField force;
  GaugeField1D tmp; 
  GaugeField1D v, w, c;
  GaugeField1D WupMu, VupNu;

  double* v_ptr     = v.data.getaddr(0);
  double* w_ptr     = w.data.getaddr(0);
  double* VupNu_ptr = VupNu.data.getaddr(0);
  double* WupMu_ptr = WupMu.data.getaddr(0);
  double* c_ptr     = c.data.getaddr(0);
  double* tmp_ptr   = tmp.data.getaddr(0);

  for(int m = 0; m < NDIM_; ++m){
    tmp = 0.0;
    for(int n=0; n< NDIM_; ++n){
      if(n != m){
	//Explicit staple calculation avoiding temporaries

	DirSlice_ALIGNED(v, *u_, m);
	DirSlice_ALIGNED(w, *u_, n);

	shiftField(VupNu,v_ptr,n,Forward());
	shiftField(WupMu,w_ptr,m,Forward());
	
	SRWilsonSU3_MatMult_NND(c_ptr, w_ptr, VupNu_ptr, WupMu_ptr, Nvol_);
	tmp += c;
	
	// temporal hack
	SRWilsonSU3_MatMult_NN(c_ptr, v_ptr, WupMu_ptr, Nvol_);
	SRWilsonSU3_MatMult_DN(VupNu_ptr, w_ptr,c_ptr , Nvol_);
	//SRWilsonSU3_MatMult_DNN(VupNu_ptr, w_ptr, v_ptr, WupMu_ptr, Nvol_);

	shiftField(w,VupNu_ptr,n,Backward());
	tmp += w;
      }
    }
    SRWilsonSU3_MatMult_ND(c_ptr,u_->data.getaddr(0)+18*Nvol_*m,tmp_ptr,Nvol_);
    SetSlice(force, TracelessAntihermite(c), m);
  }
  force *= 0.5*Params.beta/NC_;

  _MonitorMsg(ACTION_VERB_LEVEL, Action, force, "ActionGaugeWilson");
  return force;
}


