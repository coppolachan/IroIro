/*! @file test_ChiralCondensate.hpp
 *  @brief Declaration of classes for testing the ChiralCondensate classes
 */
#ifndef TEST_CHIRALCOND_INCLUDED
#define TEST_CHIRALCOND_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"
class Test_ChiralCondensate{
  const Measurements::Input input_;
 public:
  Test_ChiralCondensate(const Measurements::Input& input):input_(input){
  }
  int run();
};

#endif
