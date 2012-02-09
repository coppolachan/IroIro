//------------------------------------------------------------------------
/*!
 * @file testerDWF.cpp 
 * @brief Main source code for testing the DomainWall classes
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "test_DomainWall.hpp"

using namespace XML;

int main(){

  
  //Reading input file
  node top_node = getInputXML("test_DomainWall.xml");  

  //Initializing geometry using XML input
  Geometry geom(top_node);

  //Initialize GaugeField using XML input
  GaugeField GaugeF(geom);
  GaugeF.initialize(top_node);

  /////////////
  
  node DWF_node = top_node;
  descend(DWF_node, "TestOptimalDomainWall");
  
  
  Test_optimalDomainWall DomWallTest(DWF_node,GaugeF);
  DomWallTest.run();

  return 0;
}

