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

  /************************************************************************************/
  //
  // For 5-D Inversion
  //
  XML::node Kernel_node = gw_node;
  XML::descend(Kernel_node,"KernelDWF_4d");
  auto_ptr<DiracDWF4dFactory> 
    Wilson_Kernel_4d_factory(Diracs::createDiracDWF4dFactory(Kernel_node));
  auto_ptr<Dirac_DomainWall_4D> Wilson_Kernel_4d(Wilson_Kernel_4d_factory->getDirac(input_.config));

  
  try{
    const vector<double>& eval = input_.getEmodes()->evals;
    const vector<Field>& evec = input_.getEmodes()->evecs;

    int Neig = evec.size();
    CCIO::cout<<"CP:Neig= "<<Neig<<endl;

    double mass;
    XML::read(gw_node,"mass", mass,MANDATORY);
    Dirac_Wilson Dw(mass,input_.getGconf());
    CCIO::cout<<"CP:Pairing:mass="<< mass <<endl;

    for(int i=0; i<Neig-1; ++i){
      double sign = eval[i]/sqrt(eval[i]*eval[i]);
      double lmd = 0.0;
      if(abs(eval[i])<0.001)
	CCIO::cout<<"CP:lmd=0 beacause of near 0-mode eval["<<i<<"]="
		  <<eval[i]<<endl;
      else
	lmd = sign*sqrt(abs((eval[i]*eval[i] -mass*mass)/(1.0 -mass*mass)));

      Field phi = Wilson_Kernel_4d->gamma5(evec[i]);
      Field InvDg5 = Wilson_Kernel_4d->mult_inv(evec[i]);
   
      
      double fmij = phi*evec[i];
      fmij -= eval[i]/(1.0-mass); 
      fmij += phi*evec[i]*mass/(1.0-mass);

      double beta = mass/(eval[i]+(1.0 -mass)*lmd);

      CCIO::cout<<"CP:Pairing: "<< i <<" sign "<< sign <<endl;
      CCIO::cout<<"CP:Pairing: "<< i <<" beta "<< beta <<endl;
      CCIO::cout<<"CP:Pairing: "<< i <<" massless_eval "<< lmd <<endl;

      double fii = -fmij*(1.0 +2.0*beta*lmd+beta*beta)
	      +2.0*beta -2.0*beta*lmd*lmd +2.0*beta*beta*lmd
              -2.0*beta*beta*lmd*lmd*lmd;
      fii /= 2.0*beta*fmij -1.0 +2.0*beta*lmd -beta*beta*(1.0 -2.0*lmd*lmd);

      double g5fact = (eval[i]*eval[i]+mass)/(eval[i]*(1.0 + mass));

      double gii = phi*evec[i]- g5fact;

      double hii = phi*InvDg5 - (2*g5fact*g5fact - 1)/eval[i];


      CCIO::cout<<"CP:massive:Pairing:fmii: "<< i <<" "<< fmij <<endl;
      CCIO::cout<<"CP:massless:Pairing:fii: "<< i <<" "<< fii <<endl;
      CCIO::cout<<"CP:gii: "<< i <<" "<< gii <<endl;
      CCIO::cout<<"CP:hii: "<< i <<" "<< hii <<endl;
      CCIO::cout<<"CP:<g5 H-1 g5>: "<< i <<" "<< (phi*InvDg5) <<endl;
      CCIO::cout<<"CP:<g5>: "<< i <<" "<< (phi*evec[i]) <<endl;
      CCIO::cout<<"CP:massive_eval: "<< i << " " << eval[i] <<endl;
    }
  }catch(const char* error){
    CCIO::cout<<error<<"\n";
  }
    CCIO::cout<<"CP:Pairing:END"<<endl;

  

  return 0;
}


