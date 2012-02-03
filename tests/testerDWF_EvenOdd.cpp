//------------------------------------------------------------------------
/*!
 * @file testerDWF_EvenOdd.cpp 
 * @brief Main source code for testing the DomainWall_EvenOdd classes
 *
 * @author Jun Noaki
 */
//------------------------------------------------------------------------

#include "test_DWF_EvenOdd.hpp"

using namespace XML;

int main(){
  
  //Reading input file
  node top_node = getInputXML("test_DWF_EvenOdd.xml");  

  //Initializing geometry using XML input
  Geometry geom(top_node);

  //Initialize GaugeField using XML input
  GaugeField GaugeF(geom);
  GaugeF.initialize(top_node);
  //
  
  node DWF_node_EO = top_node;
  descend(DWF_node_EO, "TestDWF_EvenOdd");
    
  Test_DWF_EvenOdd DWFTest_eo(DWF_node_EO, GaugeF);
  DWFTest_eo.run();

  return 0;
}

