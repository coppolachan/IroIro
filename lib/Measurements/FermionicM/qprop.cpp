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
  SolverOutput monitor;
  Field sol(fsize_);
  
  for(int s=0; s<Nd_;++s){
    for(int c=0; c<Nc_;++c){
      
      CCIO::cout<<" Dirac index ="<<s<<" Color ="<<c<<std::endl;
      monitor = slv_->solve(sol,D_->mult_dag(src.mksrc(s,c)));
#if VERBOSITY > 0
      monitor.print();
#endif
      xq.push_back(sol);
    }
  }
}

void Qprop::calc(prop_t& xq,Source& src, int Nd, int Nc) const{
  
  xq.clear();
  SolverOutput monitor;
  Field sol(fsize_);
  
  int s = Nd;
  int c = Nc;
  
  CCIO::cout<<" Dirac index ="<<s<<" Color ="<<c<<std::endl;
  monitor = slv_->solve(sol,D_->mult_dag(src.mksrc(s,c)));
#if VERBOSITY > 0
  monitor.print();
#endif
  xq.push_back(sol);
  
}

void Qprop::calc(prop_t& xq,const prop_t& prp)const{

  xq.clear();
  SolverOutput monitor;
  Field sol(fsize_);

  for(int s=0; s<Nd_;++s){
    for(int c=0; c<Nc_;++c){
      CCIO::cout <<" Dirac index ="<<s<<" Color ="<<c<< std::endl; 
      monitor = slv_->solve(sol, D_->mult_dag(prp[Nc_*s+c]));
#if VERBOSITY > 0
      monitor.print();
#endif
      xq.push_back(sol);
    }
  }
}
