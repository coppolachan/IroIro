/*! 
 * @file action_Nf2.cpp
 * @brief Definition of methods of Action_Nf2 class
 Time-stamp: <2015-03-04 11:51:11 cossu>
 */
#include "action_Nf2.hpp"
#include "Tools/randNum_MP.h"
#include "Tools/fieldUtils.hpp"
#include "include/messages_macros.hpp"

// Temporary hack
//#define ANTIPERIODIC_BC

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
  
#ifdef ANTIPERIODIC_BC
  BC->apply_bc(*u_);
#endif
  
  MPrand::mp_get_gauss(ph,rand);
  phi_= D_->mult_dag(Field(ph));
#ifdef ANTIPERIODIC_BC
  BC->apply_bc(*u_);
#endif
}

double Action_Nf2::calc_H(){ 
#ifdef ANTIPERIODIC_BC
  BC->apply_bc(*u_);
#endif
  double H_nf2 = phi_*DdagD_inv(phi_);

#ifdef ANTIPERIODIC_BC
  BC->apply_bc(*u_);
#endif

  _Message(ACTION_VERB_LEVEL, "    [Action_Nf2] H = "<< H_nf2<<"\n");
  return H_nf2;
}

GaugeField Action_Nf2::md_force(){
#ifdef ANTIPERIODIC_BC
  BC->apply_bc(*u_);
#endif
  Field eta = DdagD_inv(phi_);
  GaugeField fce(D_->md_force(eta,D_->mult(eta)));

#ifdef ANTIPERIODIC_BC
  BC->apply_bc(*u_);
#endif

  // [fce] is [U*SigmaTilde] in smearing language
  if(smeared_) smart_conf_->smeared_force(fce);
  GaugeField force = FieldUtils::TracelessAntihermite(fce);

  _MonitorMsg(ACTION_VERB_LEVEL, Action, force, "Action_Nf2");
  return force;
}
