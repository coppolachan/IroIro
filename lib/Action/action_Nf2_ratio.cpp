/*!
 * @file action_Nf2_ratio.cpp
 * @brief Definition of methods of Action_Nf2_ratio class
 */
#include "action_Nf2_ratio.hpp"
#include "Tools/fieldUtils.hpp"
#include "include/messages_macros.hpp"
#include "include/timings.hpp"

////////////// Temporary hack
//#define ANTIPERIODIC_BC 


//::::::::::::::::::::::::::::::::Observer
void Action_Nf2_ratio::observer_update() {
  D1_->update_internal_state();  
  D2_->update_internal_state();  
}

void Action_Nf2_ratio::attach_smearing(SmartConf* SmearObj) {
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

Field Action_Nf2_ratio::DdagD1_inv(const Field& src){
  Field sol(fsize_);
  SolverOutput monitor = slv1_->solve(sol,src);
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
  monitor.print();
#endif
  return sol;
}

Field Action_Nf2_ratio::DdagD2_inv(const Field& src){
  Field sol(fsize_);
  SolverOutput monitor = slv2_->solve(sol,src);
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
  monitor.print();
#endif
  return sol;
}

void Action_Nf2_ratio::init(const RandNum& rand){
  long double rnd_timing;
  std::valarray<double> ph(fsize_);

  FINE_TIMING_START(rnd_timing); 

  MPrand::mp_get_gauss(ph,rand);

  FINE_TIMING_END(rnd_timing);
  _Message(TIMING_VERB_LEVEL, "[Timing] - Action_Nf2_ratio::init"
	   << " - Random numbers timing = "
	   << rnd_timing << std::endl);     

#ifdef ANTIPERIODIC_BC
  BC->apply_bc(*u_);
#endif

  phi_= D1_->mult_dag(Field(ph));
  phi_= D2_->mult(DdagD2_inv(phi_));
  
#ifdef ANTIPERIODIC_BC
  BC->apply_bc(*u_);
#endif
}

double Action_Nf2_ratio::calc_H(){
#ifdef ANTIPERIODIC_BC
  BC->apply_bc(*u_);
#endif

  Field zeta = D2_->mult_dag(phi_);//2 flavors
  double H_nf2r = zeta*DdagD1_inv(zeta);

#ifdef ANTIPERIODIC_BC
  BC->apply_bc(*u_);
#endif

  _Message(ACTION_VERB_LEVEL,"    ["<<name_<<"] H = "<<H_nf2r<<"\n");
  return H_nf2r;
}

GaugeField Action_Nf2_ratio::md_force(){
#ifdef ANTIPERIODIC_BC
  BC->apply_bc(*u_);
#endif

  Field eta = DdagD1_inv(D2_->mult_dag(phi_));
  long double timing;
  FINE_TIMING_START(timing);

  GaugeField fce(D1_->md_force(eta,D1_->mult(eta)));
  fce -= GaugeField(D2_->md_force(eta,phi_));

  FINE_TIMING_END(timing);
  _Message(TIMING_VERB_LEVEL, "[Timing] - Action_Nf2_ratio::md_force"
           << " - Force terms timing = "
           << timing << std::endl);

#ifdef ANTIPERIODIC_BC
  BC->apply_bc(*u_);
#endif

  FINE_TIMING_START(timing);

  if(smeared_) smart_conf_->smeared_force(fce);

  FINE_TIMING_END(timing);
  _Message(TIMING_VERB_LEVEL, "[Timing] - Action_Nf2_ratio::md_force"
           << " - Smeared force timing = "
           << timing << std::endl);

  FINE_TIMING_START(timing);

  GaugeField force = FieldUtils::TracelessAntihermite(fce); 

  FINE_TIMING_END(timing);
  _Message(TIMING_VERB_LEVEL, "[Timing] - Action_Nf2_ratio::md_force"
           << " - TracelessAntihermite timing = "
           << timing << std::endl);


  _MonitorMsg(ACTION_VERB_LEVEL, Action,force, name_);
  return force;
}
