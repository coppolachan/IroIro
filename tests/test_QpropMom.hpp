/*!
 * @file test_QpropMom.hpp
 * @brief Declaration of classes for calculating Sq(p)
 */
#ifndef TEST_QPROPMOM_INCLUDED
#define TEST_QPROPMOM_INCLUDED

#include "include/common_code.hpp"
#include "tests/tests.hpp"

class Test_QpropMom: public TestGeneral{
private:
  const char* test_name;
  XML::node node_;
  GaugeField conf_;
  GaugeField smeared_u_;
public:
  Test_QpropMom(XML::node node, GaugeField conf)
    :node_(node),conf_(conf){
    test_name = "QpropMom";
    XML::descend(node_,test_name);
  }

  Test_QpropMom(GaugeField& conf):conf_(conf){}

  int run();
};

#endif
