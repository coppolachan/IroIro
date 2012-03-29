//------------------------------------------------------------------------
/*!
 * @file testerHMCeven-odd.cpp 
 * @brief Main source code for testing the %HMC Even Odd version class
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "test_HMC_EvenOdd.hpp"
#include "include/commandline.hpp"

using namespace XML;

int main(int argc, char* argv[]){

  CommandOptions Options = ReadCmdLine(argc, argv);
  
  //Reading input file
  node top_node = getInputXML(Options.filename);  

  //Initializing geometry using XML input
  Geometry geom(top_node);

  //Initialize GaugeField using XML input
  GaugeGlobal GaugeF(geom);
  GaugeF.initialize(top_node);
  /////////////
  
  node HMC_node = top_node;
  descend(HMC_node, "HMC");
  
  
  Test_HMC_EvenOdd HMCTest(HMC_node,GaugeF);
  HMCTest.run();

  return 0;
}
