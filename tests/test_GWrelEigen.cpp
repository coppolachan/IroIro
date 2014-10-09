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
    const vector<double>& eval = input_.getEmodes()->evals;
    const vector<Field>& evec = input_.getEmodes()->evecs;

    int Neig = evec.size();
    CCIO::cout<<"CP:Neig= "<<Neig<<endl;

    double mass;
    XML::read(gw_node,"mass", mass,MANDATORY);
    Dirac_Wilson Dw(mass,input_.getGconf());
    CCIO::cout<<"CP:Pairling:mass="<< mass <<endl;

    for(int i=0; i<Neig; ++i){
      double sign = eval[i]/sqrt(eval[i]*eval[i]);
      double lmd = 0.0;
      if(abs(eval[i])<0.001)
	CCIO::cout<<"CP:lmd=0 beacause of near 0-mode eval["<<i<<"]="
		  <<eval[i]<<endl;
      else
	lmd = sign*sqrt(abs((eval[i]*eval[i] -mass*mass)/(1.0 -mass*mass)));

      Field phi = Dw.gamma5(evec[i]);
      double fmij = phi*evec[i];
      fmij -= eval[i]/(1.0-mass); 
      fmij += phi*evec[i]*mass/(1.0-mass);

      double beta = mass/(eval[i]+(1.0 -mass)*lmd);

      CCIO::cout<<"CP:Pairling: "<< i <<" sign "<< sign <<endl;
      CCIO::cout<<"CP:Pairling: "<< i <<" beta "<< beta <<endl;
      CCIO::cout<<"CP:Pairling: "<< i <<" massless_eval "<< lmd <<endl;

      double fii = -fmij*(1.0 +2.0*beta*lmd+beta*beta)
	      +2.0*beta -2.0*beta*lmd*lmd +2.0*beta*beta*lmd
              -2.0*beta*beta*lmd*lmd*lmd;
      fii /= 2.0*beta*fmij -1.0 +2.0*beta*lmd -beta*beta*(1.0 -2.0*lmd*lmd);

      double gii = phi*evec[i]-(eval[i]*eval[i]+mass)
	/(eval[i]*(1.0 +mass));

      CCIO::cout<<"CP:massive:Pairling:fmii: "<< i <<" "<< fmij <<endl;
      CCIO::cout<<"CP:massless:Pairling:fii: "<< i <<" "<< fii <<endl;
      CCIO::cout<<"CP:gii: "<< i <<" "<< gii <<endl;
      CCIO::cout<<"CP:massive_eval: "<< i << " " << eval[i] <<endl;
    }
  }catch(const char* error){
    CCIO::cout<<error<<"\n";
  }
    CCIO::cout<<"CP:Pairling:END"<<endl;
  return 0;
}


