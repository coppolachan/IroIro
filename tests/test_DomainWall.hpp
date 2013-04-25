/*!
 * @file test_DomainWall.hpp
 *
 * @brief Declaration of classes for testing the Dirac_optimalDomainWall class
 *
 */
#ifndef TEST_OPTIMALDOMAINWALL_INCLUDED
#define TEST_OPTIMALDOMAINWALL_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"
#include "Dirac_ops/dirac_DomainWall.hpp"

class Test_optimalDomainWall: public TestGeneral{
private:
  const char* test_name;
  XML::node DWFnode;
  GaugeField conf_;

  int mult5d_test       (const Dirac_optimalDomainWall&,const Field&,int);
  int mult5d_dag_test   (const Dirac_optimalDomainWall&,const Field&,int);
  int mult5d_gamma5_test(const Dirac_optimalDomainWall&,const Field&,int);

public:
  Test_optimalDomainWall(XML::node node,GaugeField conf):DWFnode(node),
							 conf_(conf){
    test_name = "TestOptimalDomainWall";
    XML::descend(DWFnode, test_name, MANDATORY);
  }

  int run();
};

#endif
