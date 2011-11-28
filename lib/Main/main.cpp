/*!
 * @mainpage KEK common code for %Lattice QCD simulations
 *
 */

//------------------------------------------------------------------------
/*!
 * @file main.cpp 
 * @brief Main source code
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "include/common_code.hpp"

#include "Test/tests_all.hpp"

using namespace XML;

int main(){

  
  //Reading input file
  node top_node = getInputXML("Example.xml");  

  //Initializing geometry using XML input
  Geometry geom(top_node);

  //Initialize GaugeField using XML input
  GaugeField GaugeF(geom);
  //  GaugeF.initialize(top_node);
  GaugeF.initializeTxt("conf_04040408.txt");

  /////////////
  
  node HMC_node = top_node;
  descend(HMC_node, "HMC");

  /*
  Test_Gauge GaugeTest(GaugeF);
  GaugeTest.run();
  */

  /*
  Test_EigenModes_IRL EigenTest(GaugeF.U);
  EigenTest.run();
  */

  /*
  Test_MultiShiftSolver MShift_SolverTest(GaugeF);
  MShift_SolverTest.run();
  */ 

  /*
  Test_Overlap Overlap_Tests(GaugeF);
  Overlap_Tests.run();
  */

  /*
  Test_optimalDomainWall test_odwf(GaugeF);
  test_odwf.run();
  */

  Test_Wilson_EvenOdd test_Wilson_eo(GaugeF);
  test_Wilson_eo.run();

  /*
  Test_Wilson test_Wilson(GaugeF);
  test_Wilson.run();
  */
  return 0;
}




