/*! @file test_Polyakov.hpp
 *  @brief Declaration of classes for testing the PolyakovLoop classes
 */
#ifndef TEST_POLYAKOVLOOP_INCLUDED
#define TEST_POLYAKOVLOOP_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"
class Test_PolyakovLoop{
  const Measurements::Input input_;
 public:
  Test_PolyakovLoop(const Measurements::Input& input):input_(input){
    CCIO::cout<<"Test_PolyakovLoop called"<<std::endl;
  }
  int run();
};

#endif
