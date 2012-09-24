/*!
 * @file test_staggered.hpp
 * @brief Tests for the propagators and sources
 */
#ifndef TEST_STAGGERED_INCLUDED
#define TEST_STAGGERED_INCLUDED

#include "include/common_code.hpp"
#include "tests.hpp"

class Test_staggered : public TestGeneral{
private:
  XML::node stagg_node_;
  GaugeField& conf_;
public:
  Test_staggered(const XML::node node, GaugeField& conf)
    :stagg_node_(node),
     conf_(conf){}

  ~Test_staggered(){}
  
  int run();
};

#endif
