//------------------------------------------------------------------------
/*!
 * @file testerWilson_EvenOdd.cpp 
 * @brief Main source code for testing the Wilson_EvenOdd classes
 *
 * @author Jun Noaki
 */
//------------------------------------------------------------------------

#include "test_wilson_EvenOdd.hpp"

using namespace XML;
using namespace MapsEnv;

int main(){
  
  //Reading input file
  node top_node = getInputXML("test_Wilson_EvenOdd.xml");  

  //Initializing geometry using XML input
  Geometry geom(top_node);
  initialize_mapper();

  //Initialize GaugeField using XML input
  GaugeGlobal GaugeF(geom);
  GaugeF.initialize(top_node);

  /////////////
  node Wilson_node_EO = top_node;
  descend(Wilson_node_EO, "TestWilson_EvenOdd");
    
  Test_Wilson_EvenOdd WilsonTest_eo(Wilson_node_EO, GaugeF);
  WilsonTest_eo.run();

  return 0;
}

