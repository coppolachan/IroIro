//------------------------------------------------------------------------
/*!
 * @file test_RationalApprox.cpp
 *
 * @brief run() function for Test_RationalApprox class test
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "Tools/RationalApprox/rationalapprox.hpp"
#include "test_RationalApprox.hpp"

#include "Solver/multiShiftSolver_CG.hpp"
#include "Solver/rationalSolver.hpp"
#include "include/format_F.h"
#include "include/fopr.h"
#include "include/factories.hpp"
//#include "Measurements/FermionicM/source_types.hpp"
#include "Measurements/FermionicM/qprop_MultiShift.hpp"

#include "EigenModes/findminmax.hpp"
#include "Tools/randNum_Factory.h"
#include <assert.h>
#include <vector>

#include "Communicator/communicator.h"


using namespace std;

int Test_RationalApprox::run(){
  CCIO::cout << "Starting Rational Approximation test" << std::endl;
  
  RNG_Env::RNG = RNG_Env::createRNGfactory(RA_node);

  // Test XML constructor
  RationalApprox TestXMLApprox(RA_node);

  // Test output
  // Reconstruct and test against pow
  double x_test = 0.5;
  double exponent = TestXMLApprox.exponent();
  double reference = pow(x_test, exponent);
  
  CCIO::cout << "Reference = "<< reference << "\n";
  
  
  //Reconstruct rational expansion

  double result;
  vector<double> Res = TestXMLApprox.Residuals();
  vector<double> Poles = TestXMLApprox.Poles();
  assert(Res.size() == Poles.size());

  result = TestXMLApprox.Const();
  
  for (int i = 0; i < Res.size(); ++i) {
    result += Res[i]/(x_test + Poles[i]);
  }

  CCIO::cout << "Result = "<< result << "\n";
  CCIO::cout << "Difference = "<< result-reference << "\n";
  
 
  
  // Testing the multishift solver
  CCIO::cout << "\n";
  // Definition of source 
  prop_t  xqs;

  XML::node kernel_node = RA_node;
  XML::descend(kernel_node, "Kernel");
  DiracWilsonLikeOperatorFactory* KernelF = DiracOperators::createDiracWilsonLikeOperatorFactory(kernel_node);
  DiracWilsonLike* Kernel = KernelF->getDiracOperator(&(Gfield_.data));


  // Find Max eigenvalue
  Fopr_DdagD* FoprKernel = new Fopr_DdagD(Kernel);
  findMinMax* MinMax = new findMinMax(FoprKernel,RNG_Env::RNG->getRandomNumGenerator(),
				      Kernel->fsize());

  MinMaxOut MinMaxResult = MinMax->findExtrema();
  
  TestXMLApprox.rescale(MinMaxResult.min, MinMaxResult.max);
  
  // Definition of the Solver
  int    Niter= 5000;
  double stop_cond = 1.0e-24;
  MultiShiftSolver* Solver = 
    new MultiShiftSolver_CG(FoprKernel,
                            stop_cond,
                            Niter);
  
  RationalSolver* RASolver = new RationalSolver(Solver, TestXMLApprox);

  XML::next_sibling(kernel_node,"Source");
  SourceFactory* SrcFactory 
    = Sources::createSourceFactory<SiteIndex,Format::Format_F>(kernel_node);
  Source* src = SrcFactory->getSource(Kernel->get_fermionFormat().Nvol()*
				      Kernel->get_fermionFormat().Nex());

  
  Field solution;
  Field solution2;
  RASolver->solve(solution, src->mksrc(0,0));
  
  RASolver->solve(solution2, solution);  
  
  ////////////////////////////////
  // Check answer
  // Compare with (M^dag M)
  Field reference_sol, diff_field,temp;
  temp.resize(solution.size());
  reference_sol.resize(solution.size());
  diff_field.resize(solution.size());
  
  temp = Kernel->mult(src->mksrc(0,0));
  reference_sol = Kernel->mult_dag(temp);
  
  diff_field = reference_sol;
  diff_field -= solution2;
  
  CCIO::cout << ":::::::: Check answer -- diff (norm) = "<< diff_field.norm() <<"\n";
}
