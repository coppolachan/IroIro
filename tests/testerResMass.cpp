//------------------------------------------------------------------------
/*!
 * @file testerResMass.cpp 
 * @brief Main source code for testing the Residual Mass calculation
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "test_ResidualMass.hpp"
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
  
  node Mres_node = top_node;
  descend(Mres_node, "TestResMass");
  
  
  Test_ResMass ResMassTest(GaugeF);
  ResMassTest.run(Mres_node);

  return 0;
}

