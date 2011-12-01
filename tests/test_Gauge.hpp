/*!
  @file test_Gauge.hpp
   
  @brief Definition of Test_Gauge class for testing the gauge Field

 */

#ifndef TEST_GAUGE_INCLUDED
#define TEST_GAUGE_INCLUDED

#include "include/common_code.hpp"
#include "tests.hpp"


/*!
 *@brief Class to test the gauge field methods and measurements
 */
class Test_Gauge: public TestGeneral{
 private:
  const GaugeField& d_conf;

  int shift();
  int plaquette();
  //  int evenodd();

 public:
  Test_Gauge(GaugeField& conf):d_conf(conf){}
  int run();
};

#endif
