/*!
 * @file action_Nf2_ratio.cpp
 *
 * @brief Definition of methods of Action_Nf2_ratio class
 *
 */
#include "action_Nf2_ratio.h"
#include "include/format_F.h"
#include "Communicator/comm_io.hpp"
#include "include/macros.hpp"
#include <typeinfo> 

Field Action_Nf2_ratio::DdagD1_inv(const Field& src){
  Field sol(fsize_);
  SolverOutput monitor = slv1_->solve(sol,src);
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
  monitor.print();
#endif
  return sol;
}

Field Action_Nf2_ratio::DdagD2_inv(const Field& src){
  Field sol(fsize_);
  SolverOutput monitor = slv2_->solve(sol,src);
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
  monitor.print();
#endif
  return sol;
}

void Action_Nf2_ratio::init(const RandNum& rand,const void*){
  
  CCIO::cout<<"Action_Nf2_ratio::init"<<std::endl;
  std::valarray<double> ph(fsize_);
  
  MPrand::mp_get_gauss(ph,rand,D1_->get_fermionFormat());

  #if VERBOSITY>=DEBUG_VERB_LEVEL
  double phsum= (ph*ph).sum();
  double phnorm= Communicator::instance()->reduce_sum(phsum);
  CCIO::cout<<"ph.norm="<<sqrt(phnorm)<<std::endl;
  #endif 

  phi_= D1_->mult_dag(Field(ph));

  #if VERBOSITY>=DEBUG_VERB_LEVEL
  double phisum= phi_.norm();
  double phinorm= Communicator::instance()->reduce_sum(phisum);
  for(int i=0; i<ph.size();++i) CCIO::cout<<"phi_["<<i<<"]="
  					  << phi_[i]<<std::endl;
  #endif

  phi_= D2_->mult(DdagD2_inv(phi_));
}

double Action_Nf2_ratio::calc_H(){
  Field zeta = D2_->mult_dag(phi_);
  return zeta*DdagD1_inv(zeta);
}
  
Field Action_Nf2_ratio::md_force(const void*){
  Field eta = DdagD1_inv(D2_->mult_dag(phi_));
  Field force= D1_->md_force(eta,D1_->mult(eta));
  force -= D2_->md_force(phi_,eta);

#if VERBOSITY>=ACTION_VERB_LEVEL
  monitor_force(force, "Action_Nf2_ratio");
#endif

  return force;
}

Action_Nf2_ratio::~Action_Nf2_ratio(){}
