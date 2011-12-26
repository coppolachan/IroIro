/*! 
 * @file action_Nf2.cpp
 *
 * @brief Definition of methods of Action_Nf2 class
 *
 */

#include "action_Nf2.h"
#include "include/format_F.h"

Field Action_Nf2::DdagD_inv(const Field& src){
  int Nconv;
  double diff;
  Field sol(fsize_);
  slv_->solve(sol,src,diff,Nconv);
  return sol;
}

void Action_Nf2::init(const RandNum& rand,const void*){

  CCIO::cout<<"Action_Nf2::init"<<std::endl;
  std::valarray<double> ph(fsize_);

  MPrand::mp_get_gauss(ph,rand,D_->get_fermionFormat());
  phi_= D_->mult_dag(Field(ph));
}

double Action_Nf2::calc_H(){ 
  double H_nf2 = phi_*DdagD_inv(phi_);
  CCIO::cout << "[Action_Nf2] H = "<< H_nf2 << std::endl; 
  return H_nf2;
}

Field Action_Nf2::md_force(const void*){
  Field eta = DdagD_inv(phi_);
  Field force(D_->md_force(eta,D_->mult(eta)));
  //
  double f_re= force.average_real();
  double f_im= force.average_imag();
  CCIO::cout<<"Action_Nf2: averaged MD-force = ("
	    << f_re<<","<< f_im 
	    <<")"<< std::endl;
  double f_max= force.max_element();
  CCIO::cout<<"Action_Nf2: maximum MD-force = "
	    << f_max<< std::endl;
  //  
  return force;
}
