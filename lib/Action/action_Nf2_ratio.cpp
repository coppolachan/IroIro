/*!
 * @file action_Nf2_ratio.cpp
 *
 * @brief Definition of methods of Action_Nf2_ratio class
 */
#include "action_Nf2_ratio.hpp"
#include "include/common_fields.hpp"
#include "include/messages_macros.hpp"

FermionField Action_Nf2_ratio::DdagD1_inv(const FermionField& src){
  FermionField sol;
  SolverOutput monitor = slv1_->solve(sol.data,src.data);
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
  monitor.print();
#endif
  return sol;
}

FermionField Action_Nf2_ratio::DdagD2_inv(const FermionField& src){
  FermionField sol;
  SolverOutput monitor = slv2_->solve(sol.data,src.data);
#if VERBOSITY >= SOLV_MONITOR_VERB_LEVEL
  monitor.print();
#endif
  return sol;
}

void Action_Nf2_ratio::init(const RandNum& rand){
  std::valarray<double> ph(phi_.format.size());
  
  MPrand::mp_get_gauss(ph,rand,D1_->get_gsite(),D1_->get_fermionFormat());

  #if VERBOSITY>=DEBUG_VERB_LEVEL
  double phsum= (ph*ph).sum();
  double phnorm= Communicator::instance()->reduce_sum(phsum);
  CCIO::cout<<"[Action_Nf2_ratio::init] fsize_  ="<<fsize_<<"\n";
  CCIO::cout<<"[Action_Nf2_ratio::init] ph.norm ="<<sqrt(phnorm)<<"\n";
  #endif 
  phi_.data = D1_->mult_dag(Field(ph));

  #if VERBOSITY>=DEBUG_VERB_LEVEL
  double phisum= phi_.norm();
  double phinorm= Communicator::instance()->reduce_sum(phisum);
  #endif
  phi_.data = D2_->mult(DdagD2_inv(phi_).data);
}

double Action_Nf2_ratio::calc_H(){
  FermionField zeta;
  zeta.data = D2_->mult_dag(phi_.data);//2 flavors
  double H_Ratio = zeta.data * DdagD1_inv(zeta).data;
  _Message(ACTION_VERB_LEVEL, "    [Action_Nf2_ratio] H = " << H_Ratio <<"\n");
  return H_Ratio;
}
  
GaugeField Action_Nf2_ratio::md_force(){
  FermionField eta;
  GaugeField force;
  eta = DdagD1_inv(FermionField(D2_->mult_dag(phi_.data)));
  force.data = D1_->md_force(eta.data,D1_->mult(eta.data));
  force.data -= D2_->md_force(eta.data,phi_.data);

  GaugeField force_ta = FieldUtils::TracelessAntihermite(force); 

  _MonitorMsg(ACTION_VERB_LEVEL, Action, force_ta, "Action_Nf2_ratio");
  return force_ta;
}

Action_Nf2_ratio::~Action_Nf2_ratio(){}
