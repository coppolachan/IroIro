/*
  @file StoutSmear.hpp
  @brief Defines Stout smearing class
*/
#include "stoutSmear.hpp"
#include "Communicator/comm_io.hpp"
#include "include/messages_macros.hpp"


void Smear_Stout::BaseSmear(GaugeField& C, const GaugeField& u_in)const{
  SmearBase->smear(C, u_in);}

void Smear_Stout::derivative(GaugeField& SigmaTerm, 
			     const GaugeField& iLambda,
			     const GaugeField& Gauge)const{
  SmearBase->derivative(SigmaTerm, iLambda, Gauge);
}

double Smear_Stout::func_xi0(double w) const{
  double xi0 = sin(w)/w;
  if( w < 1e-4 ) CCIO::cout << "[Smear_stout] w too small: "<< w <<"\n";
  return xi0;
}
