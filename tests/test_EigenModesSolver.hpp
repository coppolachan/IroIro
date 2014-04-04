/*! @file test_EigenModesSolver.hpp
 *  @brief Declaration of classes for testing the EigenModesSolver classes
 */
#ifndef TEST_EIGENMODESSOLVER_INCLUDED
#define TEST_EIGENMODESSOLVER_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests.hpp"

class Test_EigenModesSolver{
  const Measurements::Input input_;
 public:
  Test_EigenModesSolver(const Measurements::Input& input):input_(input){}
  int run();
};

#endif
