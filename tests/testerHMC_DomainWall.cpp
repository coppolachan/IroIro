//------------------------------------------------------------------------
/*!
 * @file testerHMC_DomainWall.cpp 
 * @brief Main source code for testing the %HMC Domain Wall classes
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "test_HMC_DomainWall.hpp"
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
  
  node HMC_node = top_node;
  descend(HMC_node, "HMC");
  
  
  Test_HMC_DomainWall HMCTest(HMC_node,GaugeF);
  HMCTest.run();

  return 0;
}
