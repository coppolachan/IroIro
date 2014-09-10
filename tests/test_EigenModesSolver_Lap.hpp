/*! @file test_EigenModesSolver_Lap.hpp
 *  @brief Declaration of classes for testing the EigenModesSolver classes
 *  for scalar operators such as Laplacian
 */
#ifndef TEST_EIGENMODESSOLVER_LAP_INCLUDED
#define TEST_EIGENMODESSOLVER_LAP_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests.hpp"

class Test_EigenModesSolver_Lap{
  const Measurements::Input input_;
 public:
  Test_EigenModesSolver_Lap(const Measurements::Input& input):input_(input){}
  int run();
};

#endif
