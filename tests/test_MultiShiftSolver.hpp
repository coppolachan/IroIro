/*!
 * @file test_MultiShiftSolver.hpp
 *
 * @brief Class to test MultiShiftSolver_CG 
 */

#ifndef TEST_SHIFTSOLVER_INCLUDED
#define TEST_SHIFTSOLVER_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"

class Test_MultiShiftSolver : public TestGeneral{
 private:
  const char* test_name;
  XML::node MS_node;
  GaugeField Gauge;
  
  int test1();
 public:
  Test_MultiShiftSolver(XML::node node, GaugeField& conf):MS_node(node),
							  Gauge(conf)
  {
    test_name = "MultiShift";
    XML::descend(MS_node, test_name);    
  }

  int run();
};

#endif
