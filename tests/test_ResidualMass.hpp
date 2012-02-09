/*!
 * @file test_ResidualMass.hpp
 *
 * @brief Declaration of classes for calculating the Residual Mass
 *
 */
#ifndef TEST_RESMASS_INCLUDED
#define TEST_RESMASS_INCLUDED

#include "include/common_code.hpp"
#include "Dirac_ops/dirac.h"
#include "tests/tests.hpp"


class Test_ResMass: public TestGeneral{
private:
  GaugeField& conf_;

  Field delta(const Dirac_optimalDomainWall_4D*,const Field&);
public:
  Test_ResMass(GaugeField& conf):conf_(conf){}
  int run(XML::node);
};

#endif
