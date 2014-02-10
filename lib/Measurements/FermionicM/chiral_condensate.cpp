/*! 
  @file chiral_condensate.cpp
  @brief Declaration of Chiral condensate measurement class ChiralCond
 */
#include "chiral_condensate.hpp"

double ChiralCondensate::calc(Source& src, const int stochastic_noise) const{
  // Calculates the chiral condensate by 
  // stochastic estimation
  // abstract definition

  double condensate = 0.0;
  for (int i = 0; i < stochastic_noise; i++){
    // creates the source vector by summing up 
    // sources in every component
    Field global_source(fsize());
    src.refresh(); // get a new source
    
    // Fill the global_source vector
    for (int s = 0; s < ND_; s++) {
      for (int c = 0; c < NC_; c++) {
	global_source += src.mksrc(s,c);
      }
    }
    
    double gs_norm = global_source.norm();
    global_source /= gs_norm;

    // Invert the source
    Field propagated = invert(global_source);

    condensate += propagated*global_source;
  }
  condensate /= ((double)stochastic_noise*CommonPrms::instance()->Nvol());

  // Collect all results
  condensate = Communicator::instance()->reduce_sum(condensate);
  return condensate;
}

///////////////////////////////////////////////////
// For the standard dirac operator
Field ChiralCondStd::invert(Field& f)const{
  Field solution(D_->fsize());
  slv_->solve(solution, f);
  return solution;
}
int ChiralCondStd::fsize()const{
  return D_->fsize();
}

///////////////////////////////////////////////////
// For the 4D domain wall operator
// that has a simpler inverter
int ChiralCondDWF::fsize()const{
  return Ddw_.fsize();
}
Field ChiralCondDWF::invert(Field& f)const{
  return Ddw_.mult_inv(f);
}

