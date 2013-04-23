/*!
 * @file action_staggered_ratio.cpp
 * @brief Definition of methods of Action_staggered_ratio class
 */
#include "action_staggered_ratio.hpp"
#include "Tools/randNum_MP.h"
#include "Tools/fieldUtils.hpp"
#include "include/messages_macros.hpp"
#include "Fields/field_expressions.hpp"

using namespace FieldExpression;

//::::::::::::::::::::::::::::::::Observer
void Action_staggered_ratio::observer_update() {
  D_->update_internal_state();  
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

Field Action_staggered_ratio::DdagD1e_inv(const Field& src){
  Field sol(fsize_);
  SolverOutput monitor = slv1e_->solve(sol,src);
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
  monitor.print();
#endif
  return sol;
}

Field Action_staggered_ratio::DdagD1o_inv(const Field& src){
  Field sol(fsize_);
  SolverOutput monitor = slv1o_->solve(sol,src);
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
  monitor.print();
#endif
  return sol;
}

Field Action_staggered_ratio::DdagD2e_inv(const Field& src){
  Field sol(fsize_);
  SolverOutput monitor = slv2e_->solve(sol,src);
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
  monitor.print();
#endif
  return sol;
}

/*  following parts are not established */

void Action_staggered_ratio::init(const RandNum& rand){

  std::valarray<double> ph(fsize_);
  MPrand::mp_get_gauss(ph,rand);

  Field xi(ph);
  xi -= mr_*D_->mult_eo(D_->mult_oe(Field(ph)));

  MPrand::mp_get_gauss(ph,rand);
  
  xi += (mr_-1.0)*D_->mult_eo(Field(ph));
  slv2e_->solve(phi_,xi);
}

double Action_staggered_ratio::calc_H(){
  Field zeta = D_->mult_oe(phi_);
  double Hr = zeta*DdagD1o_inv(zeta);
  Hr *= mr_*mr_;
  Hr += phi_*DdagD1e_inv(phi_);

  _Message(ACTION_VERB_LEVEL,"    [Action_staggered_ratio] H = "<<Hr<<"\n");
  return Hr;
}

GaugeField Action_staggered_ratio::md_force(){

  Field psi = DdagD1e_inv(phi_);
  GaugeField fce(D_->md_force(psi,D_->mult_oe(psi)));

  psi = mr_*mr_*DdagD1o_inv(D_->mult_oe(phi_));
  fce -= D_->md_force(phi_,psi);
  psi /= mr_;
  fce -= D_->md_force(D_->mult_eo(psi),psi);

  if(smeared_) smart_conf_->smeared_force(fce);
  GaugeField force = FieldUtils::TracelessAntihermite(fce); 

  _MonitorMsg(ACTION_VERB_LEVEL, Action,force,"Action_staggered_ratio");
  return force;
}


