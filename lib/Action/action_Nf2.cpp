/*! 
 * @file action_Nf2.cpp
 *
 * @brief Definition of methods of Action_Nf2 class
 *
 */

#include "action_Nf2.h"

Field Action_Nf2::DdagD_inv(const Field& src){
  int Nconv;
  double diff;
  Field sol(fsize_);
  slv_->solve(sol,src,diff,Nconv);
  return sol;
}

void Action_Nf2::init(Field&,const RandNum& rand,const void*){
  phi_.resize(fsize_);
  std::valarray<double> ph(fsize_);
  //rand.get_gauss(ph);
  
  Format::Format_F fmt(CommonPrms::instance()->Nvol());
  MPrand::mp_get(ph,rand,fmt);

  //if(Communicator::instance()->nodeid()==0) 
  //  for(int i=0;i<ph.size();++i) std::cout<<"ph["<<i<<"]="<<ph[i]<<std::endl;
  
  phi_= D_->mult_dag(Field(ph));
}

void Action_Nf2::init(Field& P,const RandNum& rand,const Field& U, const void*){
  *u_ = U;
  Action_Nf2::init(P, rand);
}


double Action_Nf2::calc_H(){ return phi_*DdagD_inv(phi_);}

Field Action_Nf2::md_force(const void*){
  Field eta = DdagD_inv(phi_);
  return D_->md_force(eta,D_->mult(eta));
}

Field Action_Nf2::md_force(const Field& U, const void*){
    *u_ = U;
    md_force();
}

  
