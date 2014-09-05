/*! 
  @file test_QnumSuscept.hpp
  @brief calculation of building blocks of quark number susceptibility for 
  WilsonDiracFiniteDensity operators
 */
#ifndef TEST_QNUMSUSCEPT_INCLUDED
#define TEST_QNUMSUSCEPT_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests.hpp"

class Test_QnumSuscept{
  const Measurements::Input input_;
public:
  Test_QnumSuscept(const Measurements::Input& input):input_(input){}
  int run();
};

#endif
