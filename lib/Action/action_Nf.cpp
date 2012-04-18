/*! 
 * @file action_Nf2.cpp
 *
 * @brief Definition of methods of Action_Nf class
 *
 * Any number of flavours
 */
#include "action_Nf.hpp"
#include "Tools/randNum_MP.h"
#include "Tools/fieldUtils.hpp"
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


void Action_Nf::init(const RandNum& rand){
  SolverOutput monitor;
  std::valarray<double> xi(fermion_size_);

  slv_->set_Approx(PseudoFermionsApprox_);

  // Loop on pseudofermions
  for(int i=0; i<Params_.n_pseudof_; ++i){
    // Generates a gaussian distributed vector <xi>   
    MPrand::mp_get_gauss(xi,rand,D_->get_gsite(),D_->get_fermionFormat());
    
    // Generates pseudofermions <phi_> = (M^dag M)^(Nf/4n) <xi> 
    // where n is the number of pseudofermion fields 
    // and Nf the number of flavors
    monitor =  slv_->solve(phi_[i], Field(xi));
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
    monitor.print();
#endif 
  }

}

double Action_Nf::calc_H(){ 
  // Calculates action for the Metropolis step
  SolverOutput monitor;
  double H_nf2 = 0.0;
  Field temp;
  temp.resize(fermion_size_);

  slv_->set_Approx(MetropolisApprox_);

  for(int i=0; i<Params_.n_pseudof_; ++i){
    monitor = slv_->solve_inv(temp, phi_[i]); // (M^dag M)^(-Nf/2n) <phi_>
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
    monitor.print();
#endif 
    H_nf2 += (phi_[i]) * temp;
  }
  _Message(ACTION_VERB_LEVEL, "    [Action_Nf] H = "<< H_nf2 <<"\n");
  return H_nf2;
}

GaugeField Action_Nf::md_force(){
  // Calculates the force
  GaugeField fce;
  std::vector<Field> eta;
  GaugeField force;
  SolverOutput monitor;

  slv_->set_Approx(MolecularDynApprox_);

  // Loop on pseudofermions
  for (int pf = 0; pf < Params_.n_pseudof_; ++pf){ 
    
    //(M^dag M)^(-Nf/2n) <phi_>
    monitor =  slv_->solve_noReconstruct(eta, phi_[pf]); 
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
    monitor.print();
#endif 
    
    for(int i = 0; i < Params_.degree_[MDStep]; ++i) {
      
      fce.data = D_->md_force(eta[i],D_->mult(eta[i]));
      fce.data *=  MolecularDynApprox_.InvResiduals()[i];
      // [fce] is [U*SigmaTilde] in smearing language
      if(smeared_) SmartField_->smeared_force(fce);
      
      force += FieldUtils::TracelessAntihermite(fce);
      
    }
  }
  
  _MonitorMsg(ACTION_VERB_LEVEL, Action, force, "Action_Nf");
  return force;
}
