/*!
  @file action_gauge_wilson.cpp
  @brief Definition of the ActionGaugeWilson class

  Time-stamp: <2014-09-24 10:05:56 cossu>
*/
#include "action_gauge_wilson.hpp"
#include "include/messages_macros.hpp"

double ActionGaugeWilson::calc_H(){
  //Number of plaquettes
  int Nplaq = CommonPrms::instance()->Lvol()*NDIM_*(NDIM_-1)/2.0;
  double plaq = stpl_.plaquette(*u_);
  double plaq_xy = stpl_.plaq_mu_nu(*u_, 0,1);
  double plaq_xz = stpl_.plaq_mu_nu(*u_, 0,2);

  double Hgauge = Params.beta*Nplaq*(1.0 -plaq);

  _Message(ACTION_VERB_LEVEL,"    [ActionGaugeWilson] H = "<<Hgauge <<"\n");
  _Message(ACTION_VERB_LEVEL,"    [ActionGaugeWilson] Plaquette      = "<<plaq<<"\n");
  _Message(ACTION_VERB_LEVEL,"    [ActionGaugeWilson] Plaquette (xy) = "<<plaq_xy<<"\n");
  _Message(ACTION_VERB_LEVEL,"    [ActionGaugeWilson] Plaquette (xz) = "<<plaq_xz<<"\n");
  
  return Hgauge;
}

