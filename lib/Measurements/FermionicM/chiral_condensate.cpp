/*! 
  @file chiral_condensate.cpp
  @brief Declaration of Chiral condensate measurement class ChiralCond
 */
#include "chiral_condensate.hpp"

double ChiralCondStd::calc(Source& src) const{
  // Calculates the chiral condensate by 
  // stochastic extimation
  
}


double ChiralCondDWF::calc(Source& src) const{
  // Calculates the chiral condensate by 
  // stochastic extimation
  int total = 100; //now hard coded 
  
  double condensate = 0.0;
  for (int i = 0; i < total; i++){
    // creates the source vector by summing up 
    // sources in every component
    Field global_source(Ddw_.fsize());
    for (int s = 0; s < Nd_; s++) {
      for (int c = 0; c < Nc_; c++) {
	global_source += src.mksrc(s,c);
      }
    }
    
    // Invert the source
    Field propagated = Ddw_.mult_inv(global_source);
    
    condensate += propagated*global_source;
    
    
  }
  
  condensate /= ((double)total*CommonPrms::instance()->Nvol());
  //source normalization missing


  return condensate;
}

