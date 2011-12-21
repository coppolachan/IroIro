/*! 
 * @file action_Clover_Nf2_EvenOdd.cpp
 *
 * @brief Definition of methods of Action_Nf2_EvenOdd class
 *
 */

#include "action_Clover_Nf2_EvenOdd.h"

Field Action_Nf2::DdagD_inv(const Field& src){
  int Nconv;
  double diff;
  Field sol(fsize_);
  slv_->solve(sol,src,diff,Nconv);
  return sol;
}

void Action_Nf2_EvenOdd::init(const RandNum& rand,const void*){
  std::valarray<double> ph(fsize_);
  MPrand::mp_get_gauss(ph,rand,D_->get_fermionFormat());

  phi_= D_->mult_dag(Field(ph));
}

double Action_Nf2_EvenOdd::calc_H(){ return phi_*DdagD_inv(phi_);}

Field Action_Nf2::md_force(const void*){
  Field eta = DdagD_inv(phi_);
  return D_->md_force(eta,D_->mult(eta));
}


  
