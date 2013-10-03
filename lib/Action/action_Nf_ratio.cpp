/*!
 * @file action_Nf_ratio.cpp
 *
 * @brief Definition of methods of Action_Nf_ratio class
 *
 * Any number of flavours
 */
#include "action_Nf_ratio.hpp"
#include "Tools/randNum_MP.h"
#include "Tools/fieldUtils.hpp"
#include "include/messages_macros.hpp"
#include "include/timings.hpp"


//::::::::::::::::::::::::::::::::Observer
void Action_Nf_ratio::observer_update(){
  D1_->update_internal_state();
  D2_->update_internal_state();
}
//::::::::::::::::::::::::::::::::
void Action_Nf_ratio::attach_smearing(SmartConf* SmearObj) {
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
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void Action_Nf_ratio::init(const RandNum& rand){
  long double rnd_timing;
  SolverOutput monitor;
  std::valarray<double> xi(fermion_size_);
  Field temp(fermion_size_);

  slv1_->set_Approx(PseudoFermionsApprox_);
  slv2_->set_Approx(PseudoFermionsApprox_);

  // Loop on pseudofermions
  for(int i=0; i<Params_.n_pseudof_; ++i){
    // Generates a gaussian distributed vector <xi>  
    FINE_TIMING_START(rnd_timing); 

    MPrand::mp_get_gauss(xi,rand);

    FINE_TIMING_END(rnd_timing);
    _Message(TIMING_VERB_LEVEL, "[Timing] - Action_Nf_ratio::init"
	     << " - Random numbers PF#"<<i<<" timing = "
	     << rnd_timing << std::endl);   
    
    // Generates pseudofermions <phi_> = (M^dag M)^(Nf/4n) <xi> 
    // where n is the number of pseudofermion fields 
    // and Nf the number of flavors
    monitor =  slv1_->solve(temp, Field(xi));
    #if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
    monitor.print();
    #endif 
    monitor =  slv2_->solve_inv(phi_[i], temp);
    #if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
    monitor.print();
    #endif 
  }
}

double Action_Nf_ratio::calc_H(){
  // Calculates action for the Metropolis step
  SolverOutput monitor;
  double H = 0.0;
  Field temp(fermion_size_), zeta(fermion_size_);

  slv1_->set_Approx(MetropolisApprox_);
  slv2_->set_Approx(MetropolisApprox_Den_);

  for(int i=0; i<Params_.n_pseudof_; ++i){
    monitor = slv2_->solve(zeta, phi_[i]); // (M2^dag M2)^(Nf/4n) <phi_>
    #if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
    monitor.print();
    #endif     
    monitor = slv1_->solve_inv(temp, zeta); // (M1^dag M1)^(-Nf/2n) <zeta>
    #if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
    monitor.print();
    #endif     
    // temp = [ (M1^dag M1)^(-Nf/2n) ][ (M2^dag M2)^(Nf/4n) ]<phi_>

    H += zeta * temp;
  }
  _Message(ACTION_VERB_LEVEL,"    ["<<name_<<"] H = "<< H <<"\n");
  
  return H;
}


GaugeField Action_Nf_ratio::md_force(){
  using namespace FieldUtils;
  // Calculates the force
  long double force_timing;
  GaugeField fce, gauge_temp;
  Field zeta1(fermion_size_), zeta3(fermion_size_);
  Field temp(fermion_size_);
  std::vector<Field> eta1;
  std::vector<Field> eta2;
  GaugeField force;
  SolverOutput monitor;

  slv1_->set_Approx(MolecularDynApprox_);
  slv2_->set_Approx(MolecularDynApprox_Den_);

  // Loop on pseudofermions
  fce = 0.0;
  for (int pf = 0; pf < Params_.n_pseudof_; ++pf){ 
    monitor =  slv2_->solve_noReconstruct_inv(eta1, phi_[pf]); 
    #if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
    monitor.print();
    #endif 
    
    FINE_TIMING_START(force_timing); 
    //////////// Reconstruct term zeta1
    zeta1 = phi_[pf];
    zeta1 *= MolecularDynApprox_Den_.Const();
    
    for(int i = 0; i < Params_.degree_[MDStep]; ++i) {
      temp = eta1[i];
      temp *= MolecularDynApprox_Den_.Residuals()[i];
      zeta1 += temp;
    }
    //// zeta1
    FINE_TIMING_END(force_timing); 
    _Message(TIMING_VERB_LEVEL, "[Timing] - Action_Nf_ratio::md_force"
	     << " - Zeta1 PF#"<<pf<<" timing = "
	     << force_timing << std::endl);   
    
    
    monitor =  slv1_->solve_noReconstruct(eta2, zeta1); 
    #if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
    monitor.print();
    #endif     
    
    long double total_force_timing = 0.0;
    for(int i = 0; i < Params_.degree_[MDStep]; ++i) {
      gauge_temp.data = D1_->md_force(eta2[i],D1_->mult(eta2[i]));
      FINE_TIMING_START(force_timing); 
      gauge_temp.data *= MolecularDynApprox_.InvResiduals()[i];
      fce += gauge_temp;
      FINE_TIMING_END(force_timing); 
      total_force_timing += force_timing;
    }
  
    _Message(TIMING_VERB_LEVEL, "[Timing] - Action_Nf_ratio::md_force"
	     << " - force1 PF#"<<pf<<" timing = "
	     << total_force_timing << std::endl);  
    
    ///////////// Reconstruct term zeta3
    FINE_TIMING_START(force_timing); 
    zeta3 = zeta1;
    zeta3 *= MolecularDynApprox_.InvConst();
    for(int i = 0; i < Params_.degree_[MDStep]; ++i) {
      temp = eta2[i];
      temp *= MolecularDynApprox_.InvResiduals()[i];
      zeta3 += temp;
    }
    
    //// zeta3
    FINE_TIMING_END(force_timing); 
    _Message(TIMING_VERB_LEVEL, "[Timing] - Action_Nf_ratio::md_force"
	     << " - zeta3 PF#"<<pf<<" timing = "
	     << force_timing << std::endl);  
    
    
    monitor =  slv2_->solve_noReconstruct_inv(eta2, zeta3); 
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
    monitor.print();
#endif 
    
    
    FINE_TIMING_START(force_timing); 
    for(int i = 0; i < Params_.degree_[MDStep]; ++i) {
      gauge_temp.data = D2_->md_force(eta1[i],D2_->mult(eta2[i]));
      gauge_temp.data += D2_->md_force(eta2[i],D2_->mult(eta1[i]));
      gauge_temp.data *= MolecularDynApprox_Den_.Residuals()[i];
      fce += gauge_temp;
    }
    FINE_TIMING_END(force_timing); 
    _Message(TIMING_VERB_LEVEL, "[Timing] - Action_Nf_ratio::md_force"
	     << " - force2 PF#"<<pf<<" timing = "
	     << force_timing << std::endl); 
    
  }
  
  if(smeared_) smart_conf_->smeared_force(fce);
  force = FieldUtils::TracelessAntihermite(fce); 
  
  _MonitorMsg(ACTION_VERB_LEVEL, Action,force, name_);
  return force;
}


