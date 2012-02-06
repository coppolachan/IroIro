//------------------------------------------------------------------------
/*!
 * @file testerResMass.cpp 
 * @brief Main source code for testing the Residual Mass calculation
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "test_ResidualMass.hpp"

using namespace XML;

int main(){

  
  //Reading input file
  node top_node = getInputXML("test_ResMass.xml");  

  //Initializing geometry using XML input
  Geometry geom(top_node);

  //Initialize GaugeField using XML input
  GaugeField GaugeF(geom);
  GaugeF.initialize(top_node);

  /////////////
  
  node Mres_node = top_node;
  descend(Mres_node, "TestResMass");
  
  
  Test_ResMass ResMassTest(Mres_node, GaugeF);
  ResMassTest.run();

  return 0;
}

