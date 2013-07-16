/*! @file test_NoMeas.hpp
    @brief Declaration of test code which does NOTHING
*/
#ifndef TEST_NOMEAS_INCLUDED
#define TEST_NOMEAS_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"

class Test_NoMeas{
  const Measurements::Input input_;
public:
  Test_NoMeas(const Measurements::Input& input):input_(input){}
  int run(){
    CCIO::cout<<"Nothing to measure.\n";
  }
};

#endif


