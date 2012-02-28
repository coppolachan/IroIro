//------------------------------------------------------------------------
/*!
 * @file testerSmear.cpp 
 * @brief Main source code for testing the Smear classes
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "test_smear.hpp"

using namespace XML;
using namespace MapsEnv;

int main(){

  
  //Reading input file
  node top_node = getInputXML("test_Smearing.xml");  

  //Initializing geometry using XML input
  Geometry geom(top_node);
  initialize_mapper();

  //Initialize GaugeField using XML input
  GaugeGlobal GaugeF(geom);
  GaugeF.initialize(top_node);
  /////////////
  
  node Smear_node = top_node;
  descend(Smear_node, "TestSmear");
    
  Test_Smear SmearTest(Smear_node, GaugeF);
  SmearTest.run();

  return 0;
}

