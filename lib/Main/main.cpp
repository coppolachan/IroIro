/*!
 * @mainpage KEK common code for %Lattice QCD simulations
 *
 * \image html keklogo-c.jpg 
 * \image html JICFUSsymbolmark170px.jpg 
 *
 * JLQCD branch of the code for lattice simulations of QCD  
 *
 * \section Features
 *
 * Current implementation:
 * - Actions (2 flavors, 2 flavors Ratio, Overlap)
 * - %Dirac operators (Wilson, Overlap, Generalized Domain Wall (4d - 5d)
 * - Linear Solvers (Conjugate Gradient Unpreconditioned, Conjigate Gradient Preconditioned, BiConjugate Gradient)
 * - Measurements (Quark Propagator [Wilson, Domain Wall], Gauge Quantities)
 * - Random Number Generators (Mersenne Twister)
 * - %XML control of program behavior
 *
 * \authors {<a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>,  Shoji Hashimoto, Jun-Ichi Noaki}
 *
 */
#include "documentation_pages.h"
//------------------------------------------------------------------------
/*!
 * @file main.cpp 
 * @brief Main source code
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------



#include "include/common_code.hpp"

#include "tests/tests_all.hpp"

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

  node HMC_node = top_node;
  descend(HMC_node, "HMC");

  Test_HMC_DomainWall test_HMC_DomainWall(GaugeF);
  test_HMC_DomainWall.run(HMC_node);
  */

  return 0;
}




