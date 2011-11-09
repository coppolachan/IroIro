/*!
 *
 * @file qprop.cpp
 *
 * @brief Definition of Qprop member functions
 *
 */

#include "qprop.h"

void Qprop::calc(prop_t& xq,Source& src) const{

  xq.clear();
  int Nconv;
  double diff;
  Field sol(fsize_);
  
  for(int s=0; s<Nd_;++s){
    for(int c=0; c<Nc_;++c){
 
      if(Communicator::instance()->nodeid()==0) 
	std::cout<<" Dirac index ="<<s<<" Color ="<<c<<std::endl;

      slv_->solve(sol,D_->mult_dag(src.mksrc(s,c)),diff,Nconv);
      
      if(Communicator::instance()->nodeid()==0) 
	std::cout<<"s="<<s<<" c="<<c
		 <<" Nconv= "<<Nconv 
		 <<" diff= "<<diff<<std::endl;     
      xq.push_back(sol);
    }
  }
}

void Qprop::calc(prop_t& xq,const prop_t& prp)const{

  xq.clear();
  int Nconv;
  double diff;
  Field sol(fsize_);

  for(int s=0; s<Nd_;++s){
    for(int c=0; c<Nc_;++c){
      slv_->solve(sol, D_->mult_dag(prp[Nc_*s+c]),
		  diff, Nconv);

      if(Communicator::instance()->nodeid()==0) 
      std::cout<<" Dirac index ="<<s<<" Color ="<<c
	       <<" Nconv ="<<Nconv <<"diff= "<<diff<<std::endl;
      xq.push_back(sol);
    }
  }
}
