//------------------------------------------------------------------------
/*!
 * @file testerClover.cpp 
 * @brief Main source code for testing the Clover classes
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------
#include "test_clover.hpp"
#include "include/commandline.hpp"

using namespace XML;

int main(int argc, char* argv[]){

  CommandOptions Options = ReadCmdLine(argc, argv);

  //Reading input file    
  node top_node = getInputXML(Options.filename);

  //Initializing geometry using XML input
  Geometry geom(top_node);

  //Initialize GaugeField using XML input
  GaugeField GaugeF(geom);
  GaugeF.initialize(top_node);

  /////////////
  node clover_node = top_node;
  descend(clover_node, "TestClover");
    
  Test_Clover CloverTest(clover_node, GaugeF);
  CloverTest.run();

  return 0;
}

