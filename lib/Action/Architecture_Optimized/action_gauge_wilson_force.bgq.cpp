/*!
  @file action_gauge_wilson_force.bgq.cpp
  @brief Specialization of the md_force method for the ActionGaugeWilson class

  This is the BGQ optimized version

  Time-stamp: <2013-04-16 16:18:33 neo>
*/
#include "action_gauge_wilson.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/fieldUtils.hpp"
#include "include/messages_macros.hpp"

GaugeField ActionGaugeWilson::md_force(){
  using namespace FieldUtils;
  using namespace SUNmatUtils;
  using namespace Mapping;

  SUNmat pl;
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
	DirSliceBGQ(v, *u_, m);
	DirSliceBGQ(w, *u_, n);
	
	shiftField(WupMu,w_ptr,m,Forward());
	shiftField(VupNu,v_ptr,n,Forward());
	
	BGWilsonSU3_MatMult_NND(c_ptr, w_ptr, VupNu_ptr, WupMu_ptr, Nvol_);
	BGWilsonSU3_MatMult_DNN(VupNu_ptr, w_ptr, v_ptr, WupMu_ptr, Nvol_);
	shiftField(w,VupNu_ptr,n,Backward());
	tmp += c;
	tmp += w;
      }
    }
    BGWilsonSU3_MatMult_ND(c_ptr,u_->data.getaddr(0)+18*Nvol_*m,tmp_ptr,Nvol_);
    SetSlice(force, TracelessAntihermite(c), m);
  }
  force *= 0.5*Params.beta/NC_;

  _MonitorMsg(ACTION_VERB_LEVEL, Action, force, "ActionGaugeWilson");

  return force;
}


