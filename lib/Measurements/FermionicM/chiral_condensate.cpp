/*! 
  @file chiral_condensate.cpp
  @brief Declaration of Chiral condensate measurement class ChiralCond

  use with z2 noise (see hep-lat/9308105)
  implements full dilution of spin and color (see e.g. hep-lat/0505023) 
 */
#include "chiral_condensate.hpp"

double ChiralCondensate::calc(Source& src, const int stochastic_noise) const{
  // Calculates the chiral condensate by 
  // stochastic estimation (does not normalize the volume factor)
  // abstract definition
  // uses full dilution of spin and color to reduce variance

  double condensate_re = 0.0;
  double condensate_im = 0.0;
  double gamma5_top = 0.0;

  for (int i = 0; i < stochastic_noise; ++i){
    double condensate = 0.0;
    double gamma5_c = 0.0;

    //full color-spin dilution
    for (int s = 0; s< ND_; ++s) { 
      for (int c = 0; c < NC_; ++c) {  		  
	Field global_source(fsize());
	src.refresh(); // get a new random source
	global_source = src.mksrc(s,c);

	// Invert the source
	Field propagated = invert(global_source);
	Field gamma5_prop = gamma5(propagated);
	
	// Calculate traces by contracting with the source
	for (int j = 0; j < fsize(); j+=2){
	  condensate    += global_source[j]*propagated[j]  + global_source[j+1]*propagated[j+1];
	  gamma5_c      += global_source[j]*gamma5_prop[j] + global_source[j+1]*gamma5_prop[j+1];
	  condensate_im += - global_source[j+1]*propagated[j] + global_source[j]*propagated[j+1];
	}

      }// end of color dilution
    }// end of spin dilution
    
    gamma5_top    += gamma5_c;// (Tr gamma5*M^-1) \propto top_charge on each configuration
    condensate_re += condensate;// (Tr M^-1)
    
    // For systematics check 
    double temp_re = Communicator::instance()->reduce_sum(condensate_re);
    double temp_im = Communicator::instance()->reduce_sum(condensate_im);
    double temp_c  = Communicator::instance()->reduce_sum(condensate);
    double temp_g5 = Communicator::instance()->reduce_sum(gamma5_c);
    double temp_5  = Communicator::instance()->reduce_sum(gamma5_top);
    CCIO::cout << "-------- Intermediate results (not normalized by the volume)\n";
    CCIO::cout << "["<<i<<"]avg_condensate    : "<< temp_re/(i+1.0) << "   "<< temp_im/(i+1.0) <<"\n";
    CCIO::cout << "["<<i<<"]condensate        : "<< temp_c  <<"\n";
    CCIO::cout << "["<<i<<"]sq_condensate     : "<< temp_re/(i+1.0)*temp_re/(i+1.0) <<"\n";
    
    CCIO::cout << "["<<i<<"]avg_g5            : "<< temp_5/(i+1.0) <<"\n";
    CCIO::cout << "["<<i<<"]g5                : "<< temp_g5 <<"\n";
    CCIO::cout << "["<<i<<"]sq_g5             : "<< temp_5/(i+1.0)*temp_5/(i+1.0) <<"\n";
    CCIO::cout << "----------------------------------\n";
  }

  condensate_re /= ((double)stochastic_noise);
  condensate_im /= ((double)stochastic_noise);
  gamma5_top    /= ((double)stochastic_noise);

  condensate_re = Communicator::instance()->reduce_sum(condensate_re);
  condensate_im = Communicator::instance()->reduce_sum(condensate_im);
  gamma5_top    = Communicator::instance()->reduce_sum(gamma5_top);

  // Collect all results
  CCIO::cout << "-------- Results (not normalized by the volume\n";
  CCIO::cout << "V*chi     : "<< condensate_re << ","<< condensate_im << "\n";
  CCIO::cout << "q_top/m   : "<< gamma5_top << "\n";
  return condensate_re;
}

double ChiralCondensate::connected_susc(Source& src, const int stochastic_noise) const{
  // Calculates the connected susceptibilities by 
  // stochastic estimation (does not normalize the volume factor)
  // abstract definition
  // uses full dilution of spin and color to reduce variance

  double connected_ps = 0.0;
  double connected_im = 0.0;
  double connected_s = 0.0;

  for (int i = 0; i < stochastic_noise; ++i){
    double conn_ps = 0.0;
    double conn_s  = 0.0;

    //full color-spin dilution
    for (int s = 0; s< ND_; ++s) { 
      for (int c = 0; c < NC_; ++c) {  		  
	Field global_source(fsize());
	src.refresh(); // get a new random source
	global_source = src.mksrc(s,c);

	// Invert the source
	/*
	Field propagated = invert(global_source);
	Field gamma5_prop = gamma5(propagated);
	Field propagated2 = invert(propagated);
	Field propagated3 = invert(gamma5_prop);
	Field gamma5_prop2 = gamma5(propagated3);
	*/	

	Field propagated = invert(global_source);
	Field propagated2 = invert(propagated);
	Field gamma5_prop = gamma5(propagated);


	// Calculate traces by contracting with the source
	for (int j = 0; j < fsize(); j+=2){
	  conn_s  += global_source[j]*propagated2[j]  + global_source[j+1]*propagated2[j+1];
	  conn_ps += gamma5_prop[j]*gamma5_prop[j] + gamma5_prop[j+1]*gamma5_prop[j+1];
	  connected_im  += - gamma5_prop[j+1]*gamma5_prop[j] + gamma5_prop[j]*gamma5_prop[j+1];

	}

      }// end of color dilution
    }// end of spin dilution
    
    connected_ps += conn_ps; // (Tr gamma5 * M^-1 * gamma5 * M^-1)   chi_pi
    connected_s  -= conn_s;  // -(Tr M^-1 * M^-1)                     chi_delta


    // For systematics check 
    double temp_ps = Communicator::instance()->reduce_sum(connected_ps);
    double temp_im = Communicator::instance()->reduce_sum(connected_im);
    double temp_s  = Communicator::instance()->reduce_sum(connected_s);
    CCIO::cout << "-------- Intermediate results (not normalized by the volume)\n";
    CCIO::cout << "["<<i<<"]conn_ps    : "<< temp_ps/(i+1.0) << "   "<< temp_im/(i+1.0) <<"\n";
    CCIO::cout << "["<<i<<"]conn_s     : "<< temp_s/(i+1.0)  <<"\n";
    CCIO::cout << "----------------------------\n";
  }

  connected_s  /= ((double)stochastic_noise);
  connected_im /= ((double)stochastic_noise);
  connected_ps /= ((double)stochastic_noise);

  connected_ps = Communicator::instance()->reduce_sum(connected_ps);
  connected_im = Communicator::instance()->reduce_sum(connected_im);
  connected_s  = Communicator::instance()->reduce_sum(connected_s);

  // Collect all results
  CCIO::cout << "-------- Results (not normalized by the volume)\n";
  CCIO::cout << "V*chi_pi               : "<< connected_ps << ","<< connected_im << "\n";
  CCIO::cout << "V*chi_delta            : "<< connected_s << "\n";
  CCIO::cout << "V*(chi_pi - chi_delta) : "<< (connected_ps - connected_s) << "\n";
  return connected_ps;
}



///////////////////////////////////////////////////
// For the standard dirac operator
Field ChiralCondStd::invert(Field& f)const{
  Field solution(D_->fsize());
  slv_->solve(solution, f);
  return solution;
}
Field ChiralCondStd::gamma5(Field& f)const{
  return D_->gamma5(f);
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
Field ChiralCondDWF::gamma5(Field& f)const{
  return Ddw_.gamma5(f);
}
/* @brief Computes \f$ \frac{1}{1-m}(D^4_{GW}(m)^{-1}_{xy} - \delta_{xy}) */
Field ChiralCondDWF::invert(Field& f)const{
  Field inv = Ddw_.mult_inv(f);
  inv  -= f;
  inv *= one_minus_m_inv;
  return inv;
}

