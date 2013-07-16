/*!
 * @file test_wilson.hpp
 * @brief Tests for the propagators and sources
 */
#ifndef TEST_WILSON_INCLUDED
#define TEST_WILSON_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests.hpp"

class Test_Wilson{
private:
  XML::node wilson_node_;
  GaugeField& conf_;
public:
  Test_Wilson(const XML::node node, GaugeField& conf)
    :wilson_node_(node),
     conf_(conf){}

  ~Test_Wilson(){}
  
  int run();
};

#endif
