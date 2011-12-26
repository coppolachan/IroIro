//----------------------------------------------------------------------------
/*!
 * @file testerHMC_DWF.cpp 
 * @brief Main source code for testing the %HMC class for Domain Wall Fermions
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//----------------------------------------------------------------------------

#include "test_HMC_DomainWall.hpp"

using namespace XML;

int main(){

  
  //Reading input file
  node top_node = getInputXML("test_HMC_DWF.xml");  

  //Initializing geometry using XML input
  Geometry geom(top_node);

  //Initialize GaugeField using XML input
  GaugeField GaugeF(geom);
  GaugeF.initialize(top_node);

  /////////////
  
  node HMC_node = top_node;
  descend(HMC_node, "HMC");
  
  Test_HMC_DomainWall test_HMC_DomainWall(GaugeF);
  test_HMC_DomainWall.run(HMC_node);
 
  return 0;
}

