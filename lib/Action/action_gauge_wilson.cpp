/*!
  @file action_gauge_wilson.cpp
  @brief Definition of the ActionGaugeWilson class

  Time-stamp: <2013-04-16 16:16:56 neo>
*/
#include "action_gauge_wilson.hpp"
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

