/*!
  @file action_gauge_wilson.cpp

  @brief Definition of the ActionGaugeWilson class
*/
#include "action_gauge_wilson.hpp"
#include "Communicator/comm_io.hpp"
#include "Tools/sunMatUtils.hpp"

double ActionGaugeWilson::calc_H(){
  double plaq = stpl_->plaquette(*u_);
  //Number of plaquettes
  int Nplaq = CommonPrms::instance()->NP()*Nvol_*Ndim_*(Ndim_-1)/2.0;

  double Hgauge = Params.beta*Nplaq*(1.0-plaq);

  _Message(ACTION_VERB_LEVEL,"    [ActionGaugeWilson] H = "<< Hgauge <<"\n");
  CCIO::cout<<"    -- Plaquette = "<< plaq <<"\n";
  
  return Hgauge;
}


Field ActionGaugeWilson::md_force(const void*){
  using namespace SUNmat_utils;
  SUNmat pl;
  GaugeField force;
  GaugeField1D tmp; 
 
  for(int m = 0; m < Ndim_; ++m){
    tmp.U = 0.0;
    for(int n=0; n< Ndim_; ++n){
      if(n != m){
	tmp.U += stpl_->upper(*u_,m,n);
	tmp.U += stpl_->lower(*u_,m,n);
      }
    }
    for(int site=0; site<Nvol_; ++site){
      pl = (u(*u_,gf_,site,m)*u_dag(tmp,site));
      force.U.set(force.Format.cslice(0,site,m),
		  anti_hermite_traceless(pl));
    }
  }

  force.U *= 0.5*Params.beta/CommonPrms::instance()->Nc();

  _MonitorMsg(ACTION_VERB_LEVEL, Action, force.U, "ActionGaugeWilson");

  return force.U;
}


