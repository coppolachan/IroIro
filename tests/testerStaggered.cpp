//------------------------------------------------------------------------
/*!
 * @file testerStaggered.cpp 
 * @brief Main source code for testing the Dirac_staggered classes
 */
//------------------------------------------------------------------------
#include "test_staggered.hpp"

using namespace XML;

int main(){
  
  //Reading input file
  node top_node = getInputXML("test_staggered.xml");  

  //Initializing geometry using XML input
  Geometry geom(top_node);

  //Initialize GaugeField using XML input
  GaugeGlobal GaugeF(geom);
  GaugeF.initialize(top_node);

  /////////////
  node stagg_node = top_node;
  descend(stagg_node, "TestStaggered");
    
  Test_staggered staggTest(stagg_node, GaugeF);
  staggTest.run();

  return 0;
}

