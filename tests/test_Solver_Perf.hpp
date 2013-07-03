/*!
 * @file test_Solver_Perf.hpp
 *
 * @brief Declaration of classes for testing the solver and kernel performances
 *
 */
#ifndef TEST_SOLVER_PERF_INCLUDED
#define TEST_SOLVER_PERF_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"

class Test_Solver_Performance{
private:
  const char* test_name;
  XML::node DWFnode;
  GaugeField conf_;

public:
  Test_Solver_Performance(XML::node node,GaugeField conf):DWFnode(node),
							  conf_(conf){
    test_name = "Test_Performance";
    XML::descend(DWFnode, test_name, MANDATORY);
  }

  int run();
};

#endif
