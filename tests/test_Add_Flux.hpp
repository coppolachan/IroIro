/*!
  @file test_Add_Flux.hpp
  @brief Definition of Test_Add_Flux class for add an abelian Flux to the  gauge Field
 */

#ifndef TEST_ADD_FLUX_INCLUDED
#define TEST_ADD_FLUX_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests.hpp"

/*!
 *@brief Class to add the abelian Flux to the gauge field 
 */

class Test_Add_Flux{
  const char* test_name;
  XML::node Gauge_node;
  GaugeField conf_;

  int plaquette();
 public:
  Test_Add_Flux(XML::node node, GaugeField& conf)
    :Gauge_node(node),conf_(conf){
    test_name = "TestAddFlux";
    XML::descend(Gauge_node, test_name);
  }
  
  int run();
};

#endif
