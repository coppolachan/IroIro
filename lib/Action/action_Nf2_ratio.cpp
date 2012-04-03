/*!
 * @file action_Nf2_ratio.cpp
 * @brief Definition of methods of Action_Nf2_ratio class
 */
#include "action_Nf2_ratio.hpp"
#include "Tools/fieldUtils.hpp"
#include "include/messages_macros.hpp"

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

void Action_Nf2_ratio::init(const RandNum& rand){
  std::valarray<double> ph(fsize_);
  MPrand::mp_get_gauss(ph,rand,D1_->get_gsite(),D1_->get_fermionFormat());
  phi_= D1_->mult_dag(Field(ph));
  phi_= D2_->mult(DdagD2_inv(phi_));
}

double Action_Nf2_ratio::calc_H(){
  Field zeta = D2_->mult_dag(phi_);//2 flavors
  double H_nf2r = zeta*DdagD1_inv(zeta);
  _Message(ACTION_VERB_LEVEL,"    [Action_Nf2_ratio] H = "<<H_nf2r<<"\n");
  return H_nf2r;
}

void Action_Nf2_ratio::observer_update(){
  D1_->update_internal_state();
  D2_->update_internal_state();
}
  
GaugeField Action_Nf2_ratio::md_force(){
  Field eta = DdagD1_inv(D2_->mult_dag(phi_));
  GaugeField fce(D1_->md_force(eta,D1_->mult(eta)));
  fce -= GaugeField(D2_->md_force(eta,phi_));

  if(smart_conf_) smart_conf_->smeared_force(fce);
  GaugeField force = FieldUtils::TracelessAntihermite(fce); 

  _MonitorMsg(ACTION_VERB_LEVEL, Action,force,"Action_Nf2_ratio");
  return force;
}


