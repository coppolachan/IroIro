/*!
 * @file DWF_residualMass.hpp
 *
 * @brief Declaration of class to calculate the Residual Mass
 * 
 *
 */
#ifndef DWF_RESMASS_HPP_
#define DWF_RESMASS_HPP_

#include "include/common_fields.hpp"
#include "Dirac_ops/dirac_DomainWall_4D.hpp"

class DWFresidualMass {
  GaugeField& conf_;
  //Currently defined only on DomainWallFermions
  Field delta(const Dirac_optimalDomainWall_4D* DWF, const Field& phi);

public:
  double calc();

};


#endif
