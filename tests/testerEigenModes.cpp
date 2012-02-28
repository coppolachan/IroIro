//------------------------------------------------------------------------
/*!
 * @file testerEigenModes.cpp 
 * @brief Main source code for testing the EigenModes classes
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "test_EigenModes_IRL.hpp"
#include "include/commandline.hpp"

using namespace XML;
using namespace MapsEnv;

int main(int argc, char* argv[]){

  CommandOptions Options = ReadCmdLine(argc, argv);
  
  //Reading input file
  node top_node = getInputXML(Options.filename);  

  //Initializing geometry using XML input
  Geometry geom(top_node);
  initialize_mapper();

  //Initialize GaugeField using XML input
  GaugeGlobal GaugeF(geom);
  GaugeF.initialize(top_node);
  /////////////
  
  node Eigen_node = top_node;
  descend(Eigen_node, "EigenTest");
  
  
  Test_EigenModes_IRL EigenTest(Eigen_node,GaugeF);
  EigenTest.run();

  return 0;
}

