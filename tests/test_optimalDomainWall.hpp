/*!
 * @file test_optimalDomainWall.hpp
 *
 * @brief Declaration of classes for testing the Dirac_optimalDomainWall class
 *
 */
#ifndef TEST_OPTIMALDOMAINWALL_INCLUDED
#define TEST_OPTIMALDOMAINWALL_INCLUDED

#include "include/common_code.hpp"
#include "Dirac_ops/dirac_optimalDomainWall.hpp"
#include "tests/tests.hpp"


class Test_optimalDomainWall: public TestGeneral{
private:
  GaugeField& conf_;

  int mult5d_test(Dirac_optimalDomainWall& DWF,
		  Field& InputField,
		  int iterations);
  int mult5d_dag_test(Dirac_optimalDomainWall& DWF,
		      Field& InputField,
		      int iterations);
  int mult5d_gamma5_test(Dirac_optimalDomainWall& DWF,
		      Field& InputField,
		      int iterations);

public:
  Test_optimalDomainWall(GaugeField& conf):conf_(conf){}
  int run(XML::node ODWFnode);
};

#endif
