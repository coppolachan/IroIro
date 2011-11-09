/*!
 * @file test_ResidualMass.hpp
 *
 * @brief Declaration of classes for calculating the Residual Mass
 *
 */
#ifndef TEST_RESMASS_INCLUDED
#define TEST_RESMASS_INCLUDED

#include "include/common_code.hpp"
#include "Dirac_ops/dirac_optimalDomainWall_4D.hpp"
#include "Test/tests.hpp"


class Test_ResMass: public TestGeneral{
private:
  GaugeField& conf_;

  const Field delta(const Dirac_optimalDomainWall_4D& DWF, Field& phi);
public:
  Test_ResMass(GaugeField& conf):conf_(conf){}
  int run(XML::node ODWFnode);
};

#endif
