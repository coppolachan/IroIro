//------------------------------------------------------------------------
/*!
 * @file testerOverlap.cpp 
 * @brief Main source code for testing the overlap classes
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "test_Overlap.hpp"
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
  
  node Overlap_node = top_node;
  descend(Overlap_node, "HMC");
  
  
  Test_Overlap OvTest(Overlap_node,GaugeF);
  OvTest.run();

  return 0;
}

