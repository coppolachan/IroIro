//------------------------------------------------------------------------
/*!
 * @file testerIO.cpp 
 * @brief Main source code for testing the input/output functions
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "test_IO.hpp"

using namespace XML;

int main(){

   
  //Reading input file
  node top_node = getInputXML("test_IO.xml");  

  //Initializing geometry using XML input
  Geometry geom(top_node);

  //Initialize GaugeField using XML input
  GaugeGlobal GaugeF(geom);
  GaugeF.initialize(top_node);

  /////////////
  
  node IO_node = top_node;
  descend(IO_node, "IOtest");
  
  
  Test_IO IOTest(IO_node,GaugeF);
  IOTest.run();
  


  return 0;
}

