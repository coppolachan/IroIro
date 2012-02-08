/*!
 * @file test_ResidualMass.hpp
 *
 * @brief Declaration of classes for calculating the Residual Mass
 *
 */
#ifndef TEST_RESMASS_INCLUDED
#define TEST_RESMASS_INCLUDED

#include "include/common_code.hpp"
#include "Dirac_ops/dirac_DomainWall_4D.hpp"
#include "tests/tests.hpp"


class Test_ResMass: public TestGeneral{
private:
  XML::node ResMassNode;
  GaugeField& conf_;
  GaugeField smeared_u_;
  GaugeField previous_u_;

  Field delta(const Dirac_optimalDomainWall_4D*,const Field&);
public:
  Test_ResMass(XML::node node, GaugeField& conf):ResMassNode(node),
						 conf_(conf){}
  int run();
};

#endif
