//------------------------------------------------------------------------
/*!
 * @file testerGauge.cpp 
 * @brief Main source code for testing the Gauge configuration classes
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "test_Gauge.hpp"
#include "include/commandline.hpp"

using namespace XML;

int main(int argc, char* argv[]){

  CommandOptions Options = ReadCmdLine(argc, argv);
  
  //Reading input file
  node top_node = getInputXML(Options.filename);  

  //Initializing geometry using XML input
  Geometry geom(top_node);

  //Initialize GaugeField using XML input
  GaugeGlobal G_Global(geom);
  G_Global.initialize(top_node);

  //GaugeField GaugeF(geom);
  //GaugeF.initialize(top_node);

  /////////////
  
  node Gauge_node = top_node;
  descend(Gauge_node, "Gauge");
  
  
  Test_Gauge GaugeTest(Gauge_node,G_Global);
  GaugeTest.run();

  return 0;
}

