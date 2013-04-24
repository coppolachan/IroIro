/*! 
 * @file action_Nf2.cpp
 * @brief Definition of methods of Action_Nf2 class
 Time-stamp: <2013-04-23 12:40:17 noaki>
 */
#include "action_Nf2.hpp"
#include "Tools/randNum_MP.h"
#include "Tools/fieldUtils.hpp"
#include "include/messages_macros.hpp"

//::::::::::::::::::::::::::::::::Observer
void Action_Nf2::observer_update() {
  D_->update_internal_state();  
}

void Action_Nf2::attach_smearing(SmartConf* SmearObj) {
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

Field Action_Nf2::DdagD_inv(const Field& src){
  Field sol(fsize_);
  SolverOutput monitor = slv_->solve(sol,src);
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
  monitor.print();
#endif
  return sol;
}

void Action_Nf2::init(const RandNum& rand){
  std::valarray<double> ph(fsize_);
  //BC->apply_bc(*u_);

  MPrand::mp_get_gauss(ph,rand);

  phi_= D_->mult_dag(Field(ph));
  //BC->apply_bc(*u_);
}

double Action_Nf2::calc_H(){ 
  //BC->apply_bc(*u_);
  double H_nf2 = phi_*DdagD_inv(phi_);
  //BC->apply_bc(*u_);
  _Message(ACTION_VERB_LEVEL, "    [Action_Nf2] H = "<< H_nf2<<"\n");
  return H_nf2;
}

GaugeField Action_Nf2::md_force(){
  //BC->apply_bc(*u_);
  Field eta = DdagD_inv(phi_);
  GaugeField fce(D_->md_force(eta,D_->mult(eta)));
  //BC->apply_bc(*u_);

  // [fce] is [U*SigmaTilde] in smearing language
  if(smeared_) smart_conf_->smeared_force(fce);
  GaugeField force = FieldUtils::TracelessAntihermite(fce);

  _MonitorMsg(ACTION_VERB_LEVEL, Action, force, "Action_Nf2");
  return force;
}
