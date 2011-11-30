//------------------------------------------------------------------------
/*!
 * @file testerODWF.cpp 
 * @brief Main source code for testing the Optimal DomainWall classes
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "test_optimalDomainWall.hpp"

using namespace XML;

int main(){

  
  //Reading input file
  node top_node = getInputXML("test_ODWF.xml");  

  //Initializing geometry using XML input
  Geometry geom(top_node);

  //Initialize GaugeField using XML input
  GaugeField GaugeF(geom);
  GaugeF.initialize(top_node);

  /////////////
  
  node ODWF_node = top_node;
  descend(ODWF_node, "TestOptimalDomainWall");
  
  
  Test_optimalDomainWall OptDomWallTest(GaugeF);
  OptDomWallTest.run(ODWF_node);

  return 0;
}

