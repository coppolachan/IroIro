/*!
 * @file test_wilson_EvenOdd.hpp
 * @brief Tests for the propagators and sources
 */
#ifndef TEST_DWF_EVENODD_INCLUDED
#define TEST_DWF_EVENODD_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests.hpp"

class Test_DWF_EvenOdd : public TestGeneral{
private:
  XML::node DWF_EO_node_;
  GaugeField& conf_;
public:

  Test_DWF_EvenOdd(const XML::node node, GaugeField& conf)
    :DWF_EO_node_(node),
     conf_(conf){}

  Test_DWF_EvenOdd(GaugeField& conf)
    :conf_(conf){}

  int run();
};

#endif
