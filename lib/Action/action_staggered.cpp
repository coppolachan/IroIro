/*! 
 * @file action_staggered.cpp
 * @brief Declaration of Action_staggered class
 */
#include "Action/action_staggered.hpp"
#include "Tools/randNum_MP.h"
#include "Tools/fieldUtils.hpp"
#include "include/messages_macros.hpp"

//::::::::::::::::::::::::::::::::Observer
void Action_staggered::observer_update(){
  D_->update_internal_state();  
}

void Action_staggered::attach_smearing(SmartConf* SmearObj) {
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

Field Action_staggered::DdagD_inv(const Field& src){
  Field sol(fsize_);
  SolverOutput monitor = slv_->solve(sol,src);
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
  monitor.print();
#endif
  return sol;
}

void Action_staggered::init(const RandNum& rand){
  std::valarray<double> ph(fsize_);
  MPrand::mp_get_gauss(ph, rand);//D_->get_RandGauss(ph,rand);
  phi_= ph;
  MPrand::mp_get_gauss(ph, rand);//???? D_->get_RandGauss(ph,rand);
  phi_-= D_->mult_eo(Field(ph));
}

double Action_staggered::calc_H(){ 
  double H = phi_*DdagD_inv(phi_);
  _Message(ACTION_VERB_LEVEL,"    [Action_staggered] H = "<<H<<"\n");
  return H;
}

GaugeField Action_staggered::md_force(){
  Field eta = DdagD_inv(phi_);
  GaugeField fce(D_->md_force(eta,D_->mult_oe(eta)));
  
  // [fce] is [U*SigmaTilde] in smearing language
  if(smeared_) smart_conf_->smeared_force(fce);
  GaugeField force = FieldUtils::TracelessAntihermite(fce);

  _MonitorMsg(ACTION_VERB_LEVEL,Action,force,"Action_staggered");
  return force;
}
