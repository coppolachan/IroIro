/*!
 * @file test_optimalDomainWall.hpp
 *
 * @brief Declaration of classes for testing the Dirac_optimalDomainWall class
 *
 */
#ifndef TEST_OPTIMALDOMAINWALL_INCLUDED
#define TEST_OPTIMALDOMAINWALL_INCLUDED

#include "include/common_code.hpp"
#include "Dirac_ops/dirac_DomainWall.hpp"
#include "tests/tests.hpp"


class Test_optimalDomainWall: public TestGeneral{
private:
  XML::node DWFnode;
  GaugeField& conf_;

  int mult5d_test(const Dirac_optimalDomainWall&,const Field&,int);
  int mult5d_dag_test(const Dirac_optimalDomainWall&,const Field&,int);
  int mult5d_gamma5_test(const Dirac_optimalDomainWall&,const Field&,int);

public:
  Test_optimalDomainWall(XML::node node,GaugeField& conf):DWFnode(node),
							  conf_(conf){}
  int run();
};

#endif
