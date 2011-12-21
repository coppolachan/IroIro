/*!
 * @file DWF_residualMass.hpp
 *
 * @brief Declaration of class to calculate the Residual Mass
 * 
 *
 */


class DWFresidualMass {
  //Currently defined only on DomainWallFermions
  const Field delta(const Dirac_DomainWall_4D* DWF, Field& phi);

public:
  double calc();

}
