/*!
 * @file test_ResidualMass.hpp
 * @brief Declaration of classes for calculating the Residual Mass
 */
#ifndef TEST_RESMASS_INCLUDED
#define TEST_RESMASS_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"
#include "Dirac_ops/dirac.hpp"

class Test_ResMass{
private:
  XML::node node_;
  GaugeField& conf_;
  GaugeField smeared_u_;
  GaugeField previous_u_;
public:
  Test_ResMass(XML::node node,GaugeField& conf)
    :node_(node),conf_(conf){}
  int run();
};

#endif
