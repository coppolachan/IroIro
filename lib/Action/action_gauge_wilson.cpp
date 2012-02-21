/*!
  @file action_gauge_wilson.cpp

  @brief Definition of the ActionGaugeWilson class
*/
#include "action_gauge_wilson.hpp"
#include "Tools/sunMatUtils.hpp"
#include "include/messages_macros.hpp"


double ActionGaugeWilson::calc_H(){
  //Number of plaquettes
  int Nplaq = CommonPrms::instance()->NP()*Nvol_*NDIM_*(NDIM_-1)/2.0;
  double plaq = stpl_.plaquette(*u_);
  double Hgauge = Params.beta * Nplaq * (1.0 - plaq);

  _Message(ACTION_VERB_LEVEL,"    [ActionGaugeWilson] H = "<< Hgauge <<"\n");
  _Message(ACTION_VERB_LEVEL,"    [ActionGaugeWilson] Plaquette = "<< plaq <<"\n");
  
  return Hgauge;
}


GaugeField ActionGaugeWilson::md_force(){
  using namespace FieldUtils;
  using namespace SUNmatUtils;
  SUNmat pl;
  GaugeField force;
  GaugeField1D tmp; 
 
  for(int m = 0; m < NDIM_; ++m){
    tmp = 0.0;
    for(int n=0; n< NDIM_; ++n){
      if(n != m){
	tmp += stpl_.upper(*u_,m,n);
	tmp += stpl_.lower(*u_,m,n);
      }
     }
     for(int site=0; site < Nvol_; ++site){
      pl = matrix(*u_, site, m) * matrix_dag(tmp, site);
      SetMatrix(force, anti_hermite(pl), site, m);
    }
  }

  force *= 0.5*Params.beta/NC_;

  _MonitorMsg(ACTION_VERB_LEVEL, Action, force, "ActionGaugeWilson");

  return force;
}


