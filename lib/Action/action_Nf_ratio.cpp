/*!
 * @file action_Nf_ratio.cpp
 * @brief Definition of methods of Action_Nf_ratio class
 */
#include "action_Nf_ratio.hpp"
#include "Tools/randNum_MP.h"
#include "Tools/fieldUtils.hpp"
#include "include/messages_macros.hpp"

Field Action_Nf_ratio::DdagD1_inv(const Field& src){
  Field sol(fermion_size_);
  SolverOutput monitor = slv1_->solve(sol,src);
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
  monitor.print();
#endif
  return sol;
}

Field Action_Nf_ratio::DdagD2_inv(const Field& src){
  Field sol(fermion_size_);
  SolverOutput monitor = slv2_->solve(sol,src);
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
  monitor.print();
#endif
  return sol;
}

void Action_Nf_ratio::init(const RandNum& rand){
  SolverOutput monitor;
  std::valarray<double> xi(fermion_size_);
  Field temp;

  slv1_->set_Approx(PseudoFermionsApprox_);
  slv2_->set_Approx(PseudoFermionsApprox_);

  // Loop on pseudofermions
  for(int i=0; i<Params_.n_pseudof_; ++i){
    MPrand::mp_get_gauss(xi,rand,D1_->get_gsite(),D1_->get_fermionFormat());

    // Generates pseudofermions <phi_> = (M^dag M)^(Nf/4n) <xi> 
    // where n is the number of pseudofermion fields 
    // and Nf the number of flavors
    monitor = slv1_->solve(temp, Field(xi)); //(D1^dag D1)^(Nf/4n)
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
    monitor.print();
#endif 
    monitor = slv2_->solve_inv(phi_[i], temp); //(D2^dag D2)^(-Nf/4n)
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
    monitor.print();
#endif 
      //phi_= D1_->mult_dag(Field(ph));
      //phi_= D2_->mult(DdagD2_inv(phi_));
  }

}

double Action_Nf_ratio::calc_H(){
  // Calculates action for the Metropolis step
  SolverOutput monitor;
  double H_nf2 = 0.0;
  Field zeta, temp;

  slv1_->set_Approx(MetropolisApprox_);
  slv2_->set_Approx(MetropolisApprox_);

  for(int i=0; i<Params_.n_pseudof_; ++i){
    monitor = slv1_->solvemult_dag(phi_);//2 flavors
    
    H_nf2r = phi[i]_ * DdagD1_inv(zeta);
    _Message(ACTION_VERB_LEVEL,"    [Action_Nf_ratio] H = "<<H_nf2r<<"\n");
  }
  return H_nf2r;
}

void Action_Nf_ratio::observer_update(){
  D1_->update_internal_state();
  D2_->update_internal_state();
}
  
GaugeField Action_Nf_ratio::md_force(){
  Field eta = DdagD1_inv(D2_->mult_dag(phi_));
  GaugeField fce(D1_->md_force(eta,D1_->mult(eta)));
  fce -= GaugeField(D2_->md_force(eta,phi_));

  if(smart_conf_) smart_conf_->smeared_force(fce);
  GaugeField force = FieldUtils::TracelessAntihermite(fce); 

  _MonitorMsg(ACTION_VERB_LEVEL, Action,force,"Action_Nf_ratio");
  return force;
}


