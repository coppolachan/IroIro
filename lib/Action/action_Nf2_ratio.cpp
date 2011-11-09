/*!
 * @file action_Nf2_ratio.cpp
 *
 * @brief Definition of methods of Action_Nf2_ratio class
 *
 */
#include "action_Nf2_ratio.h"
#include "Communicator/comm_io.hpp"

Field Action_Nf2_ratio::DdagD1_inv(const Field& src){
  int Nconv;
  double diff;
  Field sol(fsize_);
  slv1_->solve(sol,src,diff,Nconv);
  return sol;
}

Field Action_Nf2_ratio::DdagD2_inv(const Field& src){
  int Nconv;
  double diff;
  Field sol(fsize_);
  slv2_->solve(sol,src,diff,Nconv);
    return sol;
}

void Action_Nf2_ratio::init(Field&,const RandNum& rand,const void*){
  
  CCIO::cout<<"Action_ratio::init"<<std::endl;
  phi_.resize(fsize_);
  std::valarray<double> ph(fsize_);
  Format::Format_F fmt(CommonPrms::instance()->Nvol());
  MPrand::mp_get(ph,rand,fmt);
  
  double phsum= (ph*ph).sum();
  double phnorm= Communicator::instance()->reduce_sum(phsum);
  
  CCIO::cout<<"ph.norm="<<sqrt(phnorm)<<std::endl;
  
  phi_= D1_->mult_dag(Field(ph));
  phi_= D2_->mult(DdagD2_inv(phi_));
}

void Action_Nf2_ratio::init(Field& P,
			    const RandNum& rand,
			    const Field& U,
			    const void*){
  *u_ = U;
  init(P, rand);
}


double Action_Nf2_ratio::calc_H(){
  CCIO::cout<<"Action_Nf2_ratio::calc_H"<<std::endl;
  Field zeta = D2_->mult_dag(phi_);
  return zeta*DdagD1_inv(zeta);
}
  
Field Action_Nf2_ratio::md_force(const void*){
  Field eta = DdagD1_inv(D2_->mult_dag(phi_));
  Field force= D1_->md_force(eta,D1_->mult(eta));
  force -= D2_->md_force(eta,phi_);
  return force;
}

Field Action_Nf2_ratio::md_force(const Field& U, const void*){
  *u_ = U;
  md_force();
}

Action_Nf2_ratio::~Action_Nf2_ratio(){}
