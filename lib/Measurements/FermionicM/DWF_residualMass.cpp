/*!
 * @file DWF_residualMass.cpp
 * @brief Definition of class to calculate the Residual Mass
 */
#include "DWF_residualMass.hpp"

<<<<<<< HEAD
double DWFresidualMass::calc() {

  XML::descend(node, "DiracOperator");
  // operator
  // here using a specific factory since we are testing the DWF-4d operator
  DiracDWF4DfullFactory DWF_4d_Factory(node);
  Dirac_DomainWall_4D* DiracDWF_4d = DWF_4d_Factory.getDiracOperator(&(conf_.U));

  //Propagator
  QpropDWF QuarkPropagator(*DiracDWF_4d);

  vector<int> spos(4,0); 
  //Source generator
  Source_local<Format_F> src(spos,CommonPrms::instance()->Nvol());

  prop_t sq;  //Defines a vector of fields
  CCIO::cout << "Calculating propagator\n";
  QuarkPropagator.calc(sq,src);
 
  // Cycle among Dirac and color indexes and contract
  // D^-1 * Delta * D^-1
  double mres_numerator = 0;
  double im_check = 0;
  double mres_denominator = 0;
  Field Delta, Denom;
=======
Field DWFresidualMass::delta(const Dirac_optimalDomainWall_4D* DWF, const Field& phi){
  //Delta function = 1/4 * (1- sign^2(Hw))
  Field sign = DWF->signKernel(phi);    //   sign( Hw )
  Field delta = DWF->signKernel(sign);  // sign^2( Hw )
  delta -= phi;                         // sign^2( Hw ) - 1
  delta *= -0.25;                       // 1/4*(1 - sign^2( Hw ))
  return delta;
}
>>>>>>> origin/master

double DWFresidualMass::calc() {

  return 0;
}
