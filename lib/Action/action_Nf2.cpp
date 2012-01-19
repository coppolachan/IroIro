/*! 
 * @file action_Nf2.cpp
 * @brief Definition of methods of Action_Nf2 class
 */
#include "action_Nf2.hpp"
#include "include/format_F.h"

//::::::::::::::::::::::::::::::::Observer
void Action_Nf2::observer_update() {
  D_->update_internal_state();  
}
//::::::::::::::::::::::::::::::::

Field Action_Nf2::DdagD_inv(const Field& src){
  Field sol(fsize_);
  SolverOutput monitor = slv_->solve(sol,src);
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
  monitor.print();
#endif
  return sol;
}

void Action_Nf2::init(const RandNum& rand,const void*){
  std::valarray<double> ph(fsize_);
  MPrand::mp_get_gauss(ph,rand,D_->get_gsite(),D_->get_fermionFormat());
  phi_= D_->mult_dag(Field(ph));
}

double Action_Nf2::calc_H(){ 
  double H_nf2 = phi_*DdagD_inv(phi_);
  _Message(ACTION_VERB_LEVEL, "    [Action_Nf2] H = "<< H_nf2<<"\n");
  return H_nf2;
}

Field Action_Nf2::md_force(const void*){
  Field eta = DdagD_inv(phi_);
  Field force = D_->md_force(eta,D_->mult(eta));
  _MonitorMsg(ACTION_VERB_LEVEL, Action, force, "Action_Nf2");
  return force;
}
