/*!
 * @file test_Laplacian.hpp
 * @brief Test for Laplacian classes
 
 * <Time-stamp> 
 */
#ifndef TEST_LAPH_H_
#define TEST_LAPH_H_

#include "include/iroiro_code.hpp"
#include "lib/Geometry/shiftField.hpp"
class Test_LapH{
  const char* test_name;
  XML::node LapH_node;
  GaugeField Gfield_;
public:
  Test_LapH(XML::node node,GaugeField& Gfield):LapH_node(node),
					      Gfield_(Gfield){
    test_name = "LapH";
    XML::descend(LapH_node, test_name);    
  }

  int run();  
};


#endif
