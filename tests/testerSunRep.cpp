//------------------------------------------------------------------------
/*!
 * @file testerSunRep.cpp 
 * @brief Main source code for testing the sunRep class
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "test_sunRepres.hpp"

using namespace XML;

int main(){

   
  //Reading input file
  node top_node = getInputXML("test_sunRep.xml");  

  //Initializing geometry using XML input
  Geometry geom(top_node);

  //Initialize GaugeField using XML input
  GaugeGlobal GaugeF(geom);
  GaugeF.initialize(top_node);

  /////////////
  
  node Rep_node = top_node;
  descend(Rep_node, "SUNREPtest");
  
  
  Test_sunRep SUNRepTest(Rep_node,GaugeF);
  SUNRepTest.run();
  


  return 0;
}

