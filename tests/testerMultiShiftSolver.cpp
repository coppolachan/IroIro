//------------------------------------------------------------------------
/*!
 * @file testerMultishiftSolver.cpp 
 * @brief Main source code for testing the MultishiftSolver class
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "test_MultiShiftSolver.hpp"
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
  
  node MS_node = top_node;
  descend(MS_node, "MultiShift");
  
  
  Test_MultiShiftSolver MSTest(MS_node,GaugeF);
  MSTest.run();

  return 0;
}

