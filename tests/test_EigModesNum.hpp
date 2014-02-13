#ifndef TEST_EIGMODESNUM_INCLUDED
#define TEST_EIGMODESNUM_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests.hpp"

class Test_EigModesNum{
  const Measurements::Input input_;
public:
  Test_EigModesNum(const Measurements::Input& input):input_(input){}
  int run();
};

#endif
