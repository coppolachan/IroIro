/*!
 * @file test_wilson_Brillouin.hpp
 * @brief Tests for the propagators and sources
 */
#ifndef TEST_WILSON_BRILLOUIN_INCLUDED
#define TEST_WILSON_BRILLOUIN_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests.hpp"

class Test_Wilson_Brillouin{
  XML::node wilson_node_;
  GaugeField& conf_;
public:
  Test_Wilson_Brillouin(const XML::node node, GaugeField& conf)
    :wilson_node_(node),conf_(conf){}

  int run();
};

#endif
