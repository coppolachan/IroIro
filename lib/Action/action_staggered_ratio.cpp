/*!
 * @file action_staggered_ratio.cpp
 * @brief Definition of methods of Action_staggered_ratio class
 */
#include "action_staggered_ratio.hpp"
#include "Tools/randNum_MP.h"
#include "Tools/fieldUtils.hpp"
#include "include/messages_macros.hpp"

//::::::::::::::::::::::::::::::::Observer
void Action_staggered_ratio::observer_update() {
  D1_->update_internal_state();  
  D2_->update_internal_state();  
}

void Action_staggered_ratio::attach_smearing(SmartConf* SmearObj) {
  // Checks that the pointer for gauge field u_ 
  // points correctly to the smeared configuration
  // otherwise falls back to standard update
  if (u_ == SmearObj->get_current_conf()) {
    smart_conf_= SmearObj; //set the configuration
    CCIO::cout << "Succesfully attached smearing routines\n";
  }else{
    CCIO::cout << "Pointers disagree - Smearing not allowed\n";
    smeared_ = false;
  }
}

Field Action_staggered_ratio::DdagD1_inv(const Field& src){
  Field sol(fsize_);
  SolverOutput monitor = slv1_->solve(sol,src);
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
  monitor.print();
#endif
  return sol;
}

Field Action_staggered_ratio::DdagD2_inv(const Field& src){
  Field sol(fsize_);
  SolverOutput monitor = slv2_->solve(sol,src);
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
  monitor.print();
#endif
  return sol;
}

/*  following parts are not established */

void Action_staggered_ratio::init(const RandNum& rand){
  std::valarray<double> ph(fsize_);
  MPrand::mp_get_gauss(ph,rand,D1_->get_gsite(),D1_->get_fermionFormat());

  Field xi(ph);
  xi -= D2_->mult_eo(D1_->mult_oe(Field(ph)));

  MPrand::mp_get_gauss(ph,rand,D1_->get_gsite(),D1_->get_fermionFormat());

  xi += D2_->mult_eo(Field(ph));
  xi -= D1_->mult_eo(Field(ph));

  slv2_->solve(phi_,xi);
}

double Action_staggered_ratio::calc_H(){
  Field zeta = D2_->mult_dag(phi_);//2 flavors
  double H_nf2r = zeta*DdagD1_inv(zeta);
  _Message(ACTION_VERB_LEVEL,"    [Action_Nf2_ratio] H = "<<H_nf2r<<"\n");
  return H_nf2r;
}

GaugeField Action_staggered_ratio::md_force(){
  Field eta = DdagD1_inv(D2_->mult_dag(phi_));
  GaugeField fce(D1_->md_force(eta,D1_->mult(eta)));
  fce -= GaugeField(D2_->md_force(eta,phi_));

  if(smeared_) smart_conf_->smeared_force(fce);
  GaugeField force = FieldUtils::TracelessAntihermite(fce); 

  _MonitorMsg(ACTION_VERB_LEVEL, Action,force,"Action_Nf2_ratio");
  return force;
}


