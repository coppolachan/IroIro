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

/*
Field Action_Nf_ratio::DdagD1_inv(const Field& src){
  Field sol(fsize_);
  SolverOutput monitor = slv1_->solve(sol,src);
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
  monitor.print();
#endif
  return sol;
}

Field Action_Nf_ratio::DdagD2_inv(const Field& src){
  Field sol(fsize_);
  SolverOutput monitor = slv2_->solve(sol,src);
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
  monitor.print();
#endif
  return sol;
}
*/

//::::::::::::::::::::::::::::::::Observer
void Action_Nf_ratio::observer_update(){
  D1_->update_internal_state();
  D2_->update_internal_state();
}
//::::::::::::::::::::::::::::::::
void Action_Nf::attach_smearing(SmartConf* SmearObj) {}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void Action_Nf_ratio::init(const RandNum& rand){
  SolverOutput monitor;
  std::valarray<double> xi(fermion_size_);
  Field temp;
  temp.resize(fermion_size_);

  slv1_->set_Approx(PseudoFermionsApprox_);
  slv2_->set_Approx(PseudoFermionsApprox_);

  // Loop on pseudofermions
  for(int i=0; i<Params_.n_pseudof_; ++i){
    // Generates a gaussian distributed vector <xi>   
    MPrand::mp_get_gauss(xi,rand,D1_->get_gsite(),D1_->get_fermionFormat());

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
    //    phi_= D1_->mult_dag(Field(ph));
    //    phi_= D2_->mult(DdagD2_inv(phi_));
  }
}

double Action_Nf_ratio::calc_H(){
  // Calculates action for the Metropolis step
  SolverOutput monitor;
  double H_nf2 = 0.0;
  Field temp, zeta;
  temp.resize(fermion_size_);
  zeta.resize(fermion_size_);

  slv1_->set_Approx(MetropolisApprox_);
  slv2_->set_Approx(MetropolisApprox_Den_);

  for(int i=0; i<Params_.n_pseudof_; ++i){
    monitor = slv2_->solve(zeta, phi_[i]); // (M2^dag M2)^(Nf/2n) <phi_>
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
    monitor.print();
#endif     
    monitor = slv1_->solve_inv(temp, zeta); // (M1^dag M1)^(-Nf/2n) <zeta>
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
    monitor.print();
#endif     
    // temp = [ (M1^dag M1)^(-Nf/2n) ][ (M2^dag M2)^(Nf/2n) ]<phi_>

    H_nf2r += zeta * temp;
  }
  _Message(ACTION_VERB_LEVEL,"    [Action_Nf_ratio] H = "<<H_nf2r<<"\n");
  
  return H_nf2r;
}


GaugeField Action_Nf_ratio::md_force(){
  // Calculates the force
  GaugeField fce;
  Field zeta1, zeta2, zeta3;
  Field temp;
  std::vector<Field> eta1;
  std::vector<Field> eta2;
  GaugeField force;
  SolverOutput monitor;

  zeta1.resize(fermion_size_);
  zeta2.resize(fermion_size_);
  zeta3.resize(fermion_size_);
  temp.resize(fermion_size_);

  slv1_->set_Approx(MolecularDynApprox_);
  slv2_->set_Approx(MolecularDynApprox_Den_);

  // Loop on pseudofermions
  for (int pf = 0; pf < Params_.n_pseudof_; ++pf){ 
    fce.data = 0.0;

    monitor =  slv2_->solve_noReconstruct_inv(eta1, phi_[pf]); 
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
    monitor.print();
#endif 

    //Reconstruct term zeta1
    zeta1 = phi_[pf];
    zeta1 *= MolecularDynApprox_Den_.Const();
    
    for(int i = 0; i < Params_.degree_[MDStep]; ++i) {
      temp = eta1[i];
      temp *= MolecularDynApprox_Den_.Residuals()[i];
      zeta1 += temp;
    }
    //// zeta1

    monitor =  slv1_->solve_noReconstruct(eta2, zeta1); 
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
    monitor.print();
#endif     

    for(int i = 0; i < Params_.degree_[MDStep]; ++i) {
     fce.data += D1_->md_force(eta2[i],D1_->mult(eta2[i]));
    }

    //Reconstruct term zeta3
    zeta3 = temp;
    zeta3 *= MolecularDynApprox_.InvConst();
    
    for(int i = 0; i < Params_.degree_[MDStep]; ++i) {
      temp = eta2[i];
      temp *= MolecularDynApprox_.InvResiduals()[i];
      zeta3 += temp;
    }
    //// zeta3

    monitor =  slv2_->solve_noReconstruct(eta2, zeta3); 
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
    monitor.print();
#endif 


    // do we need a factor of 2????
    for(int i = 0; i < Params_.degree_[MDStep]; ++i) {
     fce.data += D2_->md_force(eta1[i],D2_->mult(eta2[i]));
    }

    if(smart_conf_) smart_conf_->smeared_force(fce);
    force += FieldUtils::TracelessAntihermite(fce); 
    
  }

  _MonitorMsg(ACTION_VERB_LEVEL, Action,force,"Action_Nf_ratio");
  return force;
}


