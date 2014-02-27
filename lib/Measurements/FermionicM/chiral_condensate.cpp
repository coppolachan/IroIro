/*! 
  @file chiral_condensate.cpp
  @brief Declaration of Chiral condensate measurement class ChiralCond

  use with z2 noise (see hep-lat/9308105)
 */
#include "chiral_condensate.hpp"

double ChiralCondensate::calc(Source& src, const int stochastic_noise) const{
  // Calculates the chiral condensate by 
  // stochastic estimation
  // abstract definition

  double condensate_re = 0.0;
  double condensate_im = 0.0;

  double condensate_sq = 0.0;
  for (int i = 0; i < stochastic_noise; i++){
    // creates the source vector by summing up 
    // sources in every component
    Field global_source(fsize());
    src.refresh(); // get a new source
    
    // Fill all components of the global_source vector
    for (int s = 0; s < ND_; s++) {
      for (int c = 0; c < NC_; c++) {
	global_source += src.mksrc(s,c);
      }
    }
    // Invert the source
    Field propagated = invert(global_source);
    
    // scalar product
    for (int j = 0; j < fsize(); j+=2){
      double temp = global_source[j]*propagated[j] + global_source[j+1]*propagated[j+1];

      condensate_sq += temp*temp;

      condensate_re += temp;
      condensate_im += - global_source[j+1]*propagated[j] + global_source[j]*propagated[j+1];
    }
    // For systematics check 
    double temp_re = Communicator::instance()->reduce_sum(condensate_re);
    double temp_im = Communicator::instance()->reduce_sum(condensate_im);
    CCIO::cout << "["<<i<<"] "<< temp_re/(i+1.0) << "   "<< temp_im/(i+1.0) <<"\n";
  }
  condensate_re = Communicator::instance()->reduce_sum(condensate_re);
  condensate_im = Communicator::instance()->reduce_sum(condensate_im);

  condensate_sq = Communicator::instance()->reduce_sum(condensate_sq);

  // Collect all results
  condensate_re /= ((double)stochastic_noise*CommonPrms::instance()->Lvol());
  condensate_im /= ((double)stochastic_noise*CommonPrms::instance()->Lvol());
  condensate_sq /= ((double)stochastic_noise*CommonPrms::instance()->Lvol());
  CCIO::cout << "c ("<< condensate_re << ","<< condensate_im << ")\n";
  CCIO::cout << "squared cond  "<< condensate_sq << "\n";
  CCIO::cout << "chi_disc: "<< condensate_sq - condensate_re*condensate_re  << "\n";

  return condensate_re;
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

/* @brief Computes \f$ \frac{1}{1-m}(D^4_{GW}(m)^{-1}_{xy} - \delta_{xy}) */
Field ChiralCondDWF::invert(Field& f)const{
  Field inv = Ddw_.mult_inv(f);
  inv  -= f;
  inv *= one_minus_m_inv;
  return inv;
}

