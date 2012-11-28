/*!
  @file action_gauge_wilson.cpp
  @brief Definition of the ActionGaugeWilson class
*/
#include "action_gauge_wilson.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/fieldUtils.hpp"
#include "include/messages_macros.hpp"

double ActionGaugeWilson::calc_H(){
  //Number of plaquettes
  int Nplaq = CommonPrms::instance()->Lvol()*NDIM_*(NDIM_-1)/2.0;
  double plaq = stpl_.plaquette(*u_);
  double Hgauge = Params.beta*Nplaq*(1.0 -plaq);

  _Message(ACTION_VERB_LEVEL,"    [ActionGaugeWilson] H = "<<Hgauge <<"\n");
  _Message(ACTION_VERB_LEVEL,"    [ActionGaugeWilson] Plaquette = "<<plaq<<"\n");
  
  return Hgauge;
}


GaugeField ActionGaugeWilson::md_force(){
  using namespace FieldUtils;
  using namespace SUNmatUtils;
  using namespace Mapping;

  SUNmat pl;
  GaugeField force;
  GaugeField1D tmp; 
 
  GaugeField1D v, w, c;
  GaugeField1D WupMu, VupNu;
  int Nvol = CommonPrms::instance()->Nvol();
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
#ifdef IBM_BGQ_WILSON 
	//Explicit staple calculation avoiding temporaries
	DirSliceBGQ(v, *u_, m);
	DirSliceBGQ(w, *u_, n);
	
	shiftField(WupMu,w_ptr ,m,Forward());
	shiftField(VupNu,v_ptr ,n,Forward());
	
	BGWilsonSU3_MatMult_NND(c_ptr , w_ptr, VupNu_ptr, WupMu_ptr, Nvol);
	BGWilsonSU3_MatMult_DNN(VupNu_ptr, w_ptr, v_ptr, WupMu_ptr, Nvol);
	shiftField(w,VupNu_ptr,n,Backward());
	tmp += c;
	tmp += w;
      }
    }
    BGWilsonSU3_MatMult_ND(c_ptr , u_->data.getaddr(0)+18*Nvol*m, tmp_ptr, Nvol);
    SetSlice(force, TracelessAntihermite(c), m);
#else
        tmp += stpl_.upper_lower(*u_,m,n);
      }
      CCIO::cout << "Tmp norm: " << tmp.norm() << " mu: " << m <<"\n";
    }


    for(int site=0; site < Nvol_; ++site){
       pl = mat(*u_,site,m)*mat_dag(tmp,site);
       SetMat(force, anti_hermite_traceless(pl), site, m);
    }
    CCIO::cout << "Force norm: " << force.norm() << " mu: " << m <<"\n";
#endif
  }

  force *= 0.5*Params.beta/NC_;

  _MonitorMsg(ACTION_VERB_LEVEL, Action, force, "ActionGaugeWilson");

  return force;
}


