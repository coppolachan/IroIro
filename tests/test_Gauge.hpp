/*!
  @file test_Gauge.hpp
  @brief Definition of Test_Gauge class for testing the gauge Field
 */

#ifndef TEST_GAUGE_INCLUDED
#define TEST_GAUGE_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests.hpp"

/*!
 *@brief Class to test the gauge field methods and measurements
 */
class Test_Gauge{
  const char* test_name;
  XML::node Gauge_node;
  GaugeField d_conf;

  int map_test();
  int plaquette();
 public:
  Test_Gauge(XML::node node, GaugeField& conf)
    :Gauge_node(node),d_conf(conf){
    test_name = "TestGauge";
    XML::descend(Gauge_node, test_name);
  }
  
  int run();
};

#endif
