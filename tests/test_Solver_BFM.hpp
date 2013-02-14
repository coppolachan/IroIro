/*!
 * @file test_Solver_BFM.hpp
 *
 * @brief Declaration of classes for testing the BFM classes
 *
 */
#ifndef TEST_SOLVER_BFM_INCLUDED
#define TEST_SOLVER_BFM_INCLUDED

#include "include/common_code.hpp"
#include "tests/tests.hpp"

class Test_Solver_BFM: public TestGeneral{
private:
  const char* test_name;
  XML::node DWFnode;
  GaugeField conf_;

public:
  Test_Solver_BFM(XML::node node,GaugeField conf):DWFnode(node),
						  conf_(conf){
    test_name = "Test_BFM";
    XML::descend(DWFnode, test_name, MANDATORY);
  }

  int run();
};

#endif
