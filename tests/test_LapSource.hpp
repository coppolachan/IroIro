#ifndef TEST_LAPSOURCE_INCLUDED
#define TEST_LAPSOURCE_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"

class Test_LapSource{
  const Measurements::Input input_;
public:
  Test_LapSource(const Measurements::Input& input):input_(input){}
  int run();
};
#endif
