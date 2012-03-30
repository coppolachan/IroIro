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

Field Action_Nf::DdagD_inv(const Field& src){
  // Calculates the vector
  // v = (M^dag M)^(-Nf/2n) x

  Field sol(fermion_size_);

  SolverOutput monitor = slv_->solve(sol,src);
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
  monitor.print();
#endif
  return sol;
}

void Action_Nf::init(const RandNum& rand){
  std::valarray<double> xi(fermion_size_);
  slv_->set_Approx(PseudoFermionsApprox_);

  // Loop on pseudofermions
  for(int i=0; i<Params_.n_pseudof_; ++i){
    // Generates a gaussian distributed vector <xi>   
    MPrand::mp_get_gauss(xi,rand,D_->get_gsite(),D_->get_fermionFormat());
    
    // Generates pseudofermions <phi_> = (M^dag M)^(Nf/4n) <xi> 
    // where n is the number of pseudofermion fields 
    // and Nf the number of flavors
    slv_->solve(phi_[i], Field(xi)); 
  }

}

double Action_Nf::calc_H(){ 
  // Calculates action for the Metropolis step
  double H_nf2 = 0.0;
  Field temp;
  slv_->set_Approx(MetropolisApprox_);
  for(int i=0; i<Params_.n_pseudof_; ++i){
    slv_->solve_inv(temp, phi_[i]); // (M^dag M)^(-Nf/2n) <phi_>
    H_nf2 += phi_[i]*temp;
  }
  _Message(ACTION_VERB_LEVEL, "    [Action_Nf] H = "<< H_nf2<<"\n");
  return H_nf2;
}

GaugeField Action_Nf::md_force(){
  // Calculates the force

  Field eta(fermion_size_), eta2(fermion_size_);
  GaugeField force =0.0;

  slv_->set_Approx(MolecularDynApprox_);

  // Loop on pseudofermions
  for (int pf=0; pf < phi_.size(); ++pf){ 

    slv_->solve_inv(eta,phi_[pf]); //(M^dag M)^(-Nf/2n) <phi_>
    slv_->solve(eta2,eta);  //(M^dag M)^(Nf/2n) <phi_>

    GaugeField fce(D_->md_force(eta,eta2));
    // [fce] is [U*SigmaTilde] in smearing language
    if(smeared_) SmartField_->smeared_force(fce);

    force += FieldUtils::TracelessAntihermite(fce);
  }

  _MonitorMsg(ACTION_VERB_LEVEL, Action, force, "Action_Nf");
  return force;
}
