/*!
 * @file DWF_residualMass.cpp
 *
 * @brief Definition of class to calculate the Residual Mass
 */
#include "DWF_residualMass.hpp"

Field DWFresidualMass::delta(const Dirac_optimalDomainWall_4D* DWF, const Field& phi){
  //Delta function = 1/4 * (1- sign^2(Hw))
  Field sign = DWF->signKernel(phi);    //   sign( Hw )
  Field delta = DWF->signKernel(sign);  // sign^2( Hw )
  delta -= phi;                         // sign^2( Hw ) - 1
  delta *= -0.25;                       // 1/4*(1 - sign^2( Hw ))
  return delta;
}

double DWFresidualMass::calc() {

  return 0;
}
