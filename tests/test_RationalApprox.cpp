//------------------------------------------------------------------------
/*!
 * @file test_RationalApprox.cpp
 * @brief run() function for Test_RationalApprox class test
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "Tools/RationalApprox/rationalapprox.hpp"
#include "test_RationalApprox.hpp"
#include "Dirac_ops/dirac_Operator_FactoryCreator.hpp"
#include "Solver/multiShiftSolver_CG.hpp"
#include "Solver/rationalSolver_CG.hpp"
#include "include/format_F.h"
#include "include/fopr.h"
#include "include/factories.hpp"
#include "Measurements/FermionicM/qprop_MultiShift.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include "Smearing/smearingFactories.hpp"

#include "EigenModes/findminmax.hpp"
#include "Tools/randNum_Factory.hpp"
#include <assert.h>
#include <vector>

using namespace std;

int Test_RationalApprox::run(){
  CCIO::cout << "Starting Rational Approximation test "
	     <<"with smearing support (devel)\n";
  
  //RNG_Env::RNG = RNG_Env::createRNGfactory(RA_node);
  RNG_Env::initialize(RA_node);

// Prints plaquette (thin) link
  Staples Staple;
  CCIO::cout << "Plaquette (thin): " << Staple.plaquette(Gfield_) << std::endl;
  
///////////////////////////////////////////////////////////////////////////////
  // Smearing objects
  GaugeField smeared_u_;
  GaugeField previous_u_;

  Smear* SmearingObj;         // Empty pointer
  int Nsmear;                 // Number of smearing steps
  XML::node SmearObjNode = RA_node;// Copy the node ("descend" function updates it)
                                 // and we want to use again for RationalApproximation
  CCIO::cout<<"XML read...\n";
  XML::descend(SmearObjNode, "Smearing");//SmearObjNode now points to <Smearing> node
  XML::read(SmearObjNode, "Nsmear", Nsmear, MANDATORY);  // Reads in <Nsmear>
  
  // Create smearing factory from node information
  CCIO::cout<<"SmearingFactory...\n";
  SmearingFactory* Sm_Factory = 
    Smearings::createSmearingFactory(SmearObjNode);
  
  CCIO::cout<<"SmearingObject...\n";
  // Create smearing objects from the factory
  SmearingObj = Sm_Factory->getSmearing();
  
  // Copy original configuration to smeared_u_ 
  // smeared_u_ will be passed to the operators
  smeared_u_ = Gfield_;

  
  // Do the actual smearing 
  for(int i=0; i<Nsmear; i++) {
    previous_u_= smeared_u_;
    SmearingObj->smear(smeared_u_,previous_u_);
  }
  CCIO::cout<<"Plaquette (smeared): "<< Staple.plaquette(smeared_u_)<<std::endl;
  //////////////////////////////////////////////////////////////////////////////
  

  // Test XML constructor
  RationalApprox TestXMLApprox(RA_node);

  // Test output
  CCIO::cout << " -- Simple numerical test - check against C functions\n";
  // Reconstruct and test against pow
  double x_test = 0.5;
  double exponent = TestXMLApprox.exponent();
  double reference = pow(x_test, exponent);

  CCIO::cout << "Reference pow("<<x_test<<","<<exponent<<") = "
	     << reference << "\n";
   
  //Reconstruct rational expansion
  double result;
  vector<double> Res = TestXMLApprox.Residuals();
  vector<double> Poles = TestXMLApprox.Poles();
  assert(Res.size() == Poles.size());

  result = TestXMLApprox.Const();
  
  for (int i = 0; i < Res.size(); ++i) {
    result += Res[i]/(x_test + Poles[i]);
  }

  CCIO::cout << "Rational Approximation result = "<< result << "\n";
  CCIO::cout << "Difference = "<< result-reference << "\n";
  
  // Testing the multishift solver
  CCIO::cout << "------------------------ \n";
  // Definition of source 
  prop_t  xqs;

  XML::node kernel_node = RA_node;
  XML::descend(kernel_node, "Kernel");
  DiracWilsonLikeFactory* KernelF = 
    Diracs::createDiracWilsonLikeFactory(kernel_node);
  InputConfig input(&smeared_u_);
  DiracWilsonLike* Kernel = KernelF->getDirac(input);

  // Find Max eigenvalue
  Fopr_DdagD* FoprKernel = new Fopr_DdagD(Kernel);
  findMinMax* MinMax = new findMinMax(FoprKernel,
				      RNG_Env::RandNumG::instance().getRNG(),
				      Kernel->fsize());

  MinMaxOut MinMaxResult = MinMax->findExtrema();
  
  CCIO::cout << "Rescaling rational approximation range... ";
  TestXMLApprox.rescale(MinMaxResult.min, MinMaxResult.max);
  CCIO::cout << "done\n";

  // Definition of the Solver
  int    Niter= 5000;
  double stop_cond = 1.0e-24;
  MultiShiftSolver* Solver = 
    new MultiShiftSolver_CG(FoprKernel,
                            stop_cond,
                            Niter);
  
  RationalSolver_CG* RASolver = new RationalSolver_CG(Solver, TestXMLApprox);

  XML::next_sibling(kernel_node,"Source");
  SourceFactory* SrcFactory 
    = Sources::createSourceFactory<SiteIndex,Format::Format_F>(kernel_node);
  Source* src = SrcFactory->getSource(Kernel->fsize()/Format::Format_F::Nin());
  
  Field solution;
  Field solution2;
  CCIO::cout << "Applying rational approximation twice... ";
  RASolver->solve(solution, src->mksrc(0,0));
  RASolver->solve(solution2, solution);  
  CCIO::cout << "done"<<std::endl;
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
  
  CCIO::cout << ":::::::: Compare with (M^dag M) - meaningful only if exponent is 1/2 \n";
  CCIO::cout << ":::::::: Check answer -- diff (norm) = "<< diff_field.norm() <<"\n";
}
