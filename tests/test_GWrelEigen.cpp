/*! @file test_GWrelEigen.cpp
 * @brief implementation of the Test_GWrelEigen class
 */
#include "test_GWrelEigen.hpp"
#include "Fields/field_expressions.hpp"
#include "Dirac_ops/dirac_wilson.hpp"

using namespace std;
using namespace FieldExpression;

int Test_GWrelEigen::run(){
  
  XML::node gw_node = input_.node;
  XML::descend(gw_node,"GWrelEigen");

  try{
    InputConfig input = input_.getConfig();
    evecs_= input_.eigen->evecs_;
    evals_= input_.eigen->evals_;

    int Neig = evecs_.size();
    CCIO::cout<<"CP:Neig= "<<Neig<<endl;

    double mass;
    XML::read(gw_node,"mass", mass,MANDATORY);
    Dirac_Wilson Dw(mass,&(input_.gconf->data));
    CCIO::cout<<"CP:Pairling:mass="<< mass <<endl;

    for(int i=0; i<Neig; ++i){
      double sign = evals_[i]/sqrt(evals_[i]*evals_[i]);
      double lmd = 0.0;
      if(abs(evals_[i])<0.001)
	CCIO::cout<<"CP:lmd=0 beacause of near 0-mode evals["<<i<<"]="
		  <<evals_[i]<<endl;
      else
	lmd = sign*sqrt(abs((evals_[i]*evals_[i] -mass*mass)/(1.0 -mass*mass)));

      Field phi = Dw.gamma5(evecs_[i]);
      double fmij = phi*evecs_[i];
      fmij -= evals_[i]/(1.0-mass); 
      fmij += phi*evecs_[i]*mass/(1.0-mass);

      double beta = mass/(evals_[i]+(1.0 -mass)*lmd);

      CCIO::cout<<"CP:Pairling: "<< i <<" sign "<< sign <<endl;
      CCIO::cout<<"CP:Pairling: "<< i <<" beta "<< beta <<endl;
      CCIO::cout<<"CP:Pairling: "<< i <<" massless_eval "<< lmd <<endl;

      double fii = -fmij*(1.0 +2.0*beta*lmd+beta*beta)
	      +2.0*beta -2.0*beta*lmd*lmd +2.0*beta*beta*lmd
              -2.0*beta*beta*lmd*lmd*lmd;
      fii /= 2.0*beta*fmij -1.0 +2.0*beta*lmd -beta*beta*(1.0 -2.0*lmd*lmd);

      double gii = phi*evecs_[i]-(evals_[i]*evals_[i]+mass)
	/(evals_[i]*(1.0 +mass));

      CCIO::cout<<"CP:massive:Pairling:fmii: "<< i <<" "<< fmij <<endl;
      CCIO::cout<<"CP:massless:Pairling:fii: "<< i <<" "<< fii <<endl;
      CCIO::cout<<"CP:gii: "<< i <<" "<< gii <<endl;
      CCIO::cout<<"CP:massive_eval: "<< i << " " << evals_[i] <<endl;
    }
  }catch(const char* error){
    CCIO::cout<<error<<"\n";
  }
    CCIO::cout<<"CP:Pairling:END"<<endl;
  return 0;
}


