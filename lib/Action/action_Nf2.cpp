/*! 
 * @file action_Nf2.cpp
 *
 * @brief Definition of methods of Action_Nf2 class
 */
#include "action_Nf2.hpp"
#include "Tools/randNum_MP.h"
#include "include/messages_macros.hpp"

//::::::::::::::::::::::::::::::::Observer
void Action_Nf2::observer_update() {
  D_->update_internal_state();  
}
//::::::::::::::::::::::::::::::::
void Action_Nf2::attach_smearing(SmartConf* SmearObj) {
  // Checks that the pointer for gauge field u_ 
  // points correctly to the smeared configuration
  // otherwise falls back to standard update
  if (u_ == SmearObj->get_current_conf()) {
    SmartField_= SmearObj; //set the configuration
    CCIO::cout << "Succesfully attached smearing routines\n";
  }else{
    CCIO::cout << "Pointers disagree - Smearing not allowed\n";
    smeared_ = false;
  }
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

FermionField Action_Nf2::DdagD_inv(const FermionField& src){
  FermionField sol;
  SolverOutput monitor = slv_->solve(sol.data,src.data);
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
  monitor.print();
#endif
  return sol;
}

void Action_Nf2::init(const RandNum& rand){
  std::valarray<double> ph(phi_.data.size());
  MPrand::mp_get_gauss(ph,rand,D_->get_gsite(),D_->get_fermionFormat());
  phi_.data = D_->mult_dag(Field(ph));
}

double Action_Nf2::calc_H(){ 
  double H_nf2 = (phi_.data) * DdagD_inv(phi_).data;
  _Message(ACTION_VERB_LEVEL, "    [Action_Nf2] H = "<< H_nf2<<"\n");
  return H_nf2;
}

GaugeField Action_Nf2::md_force(){
  FermionField eta = DdagD_inv(phi_);
  GaugeField fce;

  CCIO::cout << "ETA norm: "<< eta.norm() << std::endl;

  fce.data = D_->md_force(eta.data,D_->mult(eta.data));
  // [fce] is [U*SigmaTilde] in smearing language
  CCIO::cout << "forceD norm: "<< fce.norm() << std::endl;
  if(smeared_) SmartField_->smeared_force(fce);

  GaugeField force = FieldUtils::TracelessAntihermite(fce);

  _MonitorMsg(ACTION_VERB_LEVEL, Action, force, "Action_Nf2");
  return force;
}
