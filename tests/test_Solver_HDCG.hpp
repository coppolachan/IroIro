/*!
 * @file test_Solver_HDCG.hpp
 *
 * @brief Declaration of classes for testing the BFM HDCG classes
 *
 */
#ifndef TEST_SOLVER_BFM_HDCG_INCLUDED
#define TEST_SOLVER_BFM_HDCG_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"

class Test_Solver_HDCG: public TestGeneral{
private:
  const char* test_name;
  XML::node DWFnode;
  GaugeField conf_;

public:
  Test_Solver_HDCG(XML::node node,GaugeField conf):DWFnode(node),
						  conf_(conf){
    test_name = "Test_HDCG";
    XML::descend(DWFnode, test_name, MANDATORY);
  }

  int run();
};

#endif
