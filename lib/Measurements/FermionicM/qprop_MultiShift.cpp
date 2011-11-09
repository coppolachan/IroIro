/*!
 * @file qprop_MultiShift.cpp
 *
 * @brief Quark Propagator Qprop_MultiShift functions declaration
 *
 */

#include "qprop_MultiShift.hpp"

void Qprop_MultiShift::calc(prop_t& xq,
			    Source& Source,
			    std::vector<double> mass_shifts)const
{
 
  xq.clear();
  int Nconv;
  double diff;
  prop_t sol(Nc_*Nd_);
  for(int i=0;i<sol.size();++i) sol[i].resize(fsize_);

  for(int s=0; s<Nd_;++s){
    for(int c=0; c<Nc_;++c){
      slv_->solve(sol,Source.mksrc(s,c),mass_shifts,diff,Nconv);
      std::cout<<"Dirac Index =" <<s<< " Color ="<< c
	       <<" Iterations = "<<Nconv 
	       <<" Residual = "<<diff<<std::endl;

      for(int m=0; m<sol.size();++m) xq.push_back(D_->mult(sol[m]));
    }
  }
}
