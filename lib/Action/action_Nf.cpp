/*! 
 * @file action_Nf2.cpp
 *
 * @brief Definition of methods of Action_Nf class
 *
 * Any number of flavours
 */
#include "action_Nf.hpp"
#include "Tools/randNum_MP.h"
#include "include/messages_macros.hpp"

//::::::::::::::::::::::::::::::::Observer
void Action_Nf::observer_update() {
  D_->update_internal_state();  
}
//::::::::::::::::::::::::::::::::
void Action_Nf::attach_smearing(SmartConf* SmearObj) {
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

FermionField Action_Nf::DdagD_inv(const FermionField& src){
  FermionField sol;
  sol.resize(fermion_size_);
  SolverOutput monitor = slv_->solve(sol.data,src.data);
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
  monitor.print();
#endif
  return sol;
}

void Action_Nf::init(const RandNum& rand){
  // Generates a gaussian distributed vector <xi>
  std::valarray<double> xi(fermion_size_);
  MPrand::mp_get_gauss(xi,rand,D_->get_gsite(),D_->get_fermionFormat());
 
  // Generates pseudofermions <phi_> = (M^dag M)^(Nf/4n) <xi> 
  // where n is the number of pseudofermion fields 
  // and Nf the number of flavors
  phi_.data = D_->mult_dag(Field(ph)); 
}

double Action_Nf::calc_H(){ 
  // Calculates action for the Metropolis step
  double H_nf2 = (phi_.data) * DdagD_inv(phi_).data;
  _Message(ACTION_VERB_LEVEL, "    [Action_Nf] H = "<< H_nf2<<"\n");
  return H_nf2;
}

GaugeField Action_Nf::md_force(){
  // Calculates the force
  GaugeField fce;
  FermionField eta;
  eta.resize(fermion_size_);

  eta = DdagD_inv(phi_);

  fce.data = D_->md_force(eta.data,D_->mult(eta.data));
  // [fce] is [U*SigmaTilde] in smearing language
  if(smeared_) SmartField_->smeared_force(fce);

  GaugeField force = FieldUtils::TracelessAntihermite(fce);

  _MonitorMsg(ACTION_VERB_LEVEL, Action, force, "Action_Nf");
  return force;
}
