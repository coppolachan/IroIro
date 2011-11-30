/*!
 * @file test_MultiShiftSolver.h
 *
 * @brief Class to test MultiShiftSolver_CG 
 */

#ifndef TEST_SHIFTSOLVER_INCLUDED
#define TEST_SHIFTSOLVER_INCLUDED

#include "include/common_code.hpp"

class Test_MultiShiftSolver{
 private:
  GaugeField& Gauge;
  
  int test1();
 public:
  Test_MultiShiftSolver(GaugeField& conf):Gauge(conf){}
  int run();
};

#endif
