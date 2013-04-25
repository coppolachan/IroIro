/*!
 * @file test_wilson_Brillouin_Improved.hpp
 * @brief Tests for the propagators and sources
 */
#ifndef TEST_WILSON_BRILLOUIN_IMP_INCLUDED
#define TEST_WILSON_BRILLOUIN_IMP_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests.hpp"

class Test_Wilson_Brillouin_Imp : public TestGeneral{
private:
  XML::node wilson_node_;
  GaugeField& conf_;
public:
  Test_Wilson_Brillouin_Imp(const XML::node node, GaugeField& conf)
    :wilson_node_(node),
     conf_(conf){}

  ~Test_Wilson_Brillouin_Imp(){}
  
  int run();
};

#endif
