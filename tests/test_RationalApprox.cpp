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
#include "include/format_F.h"
#include "include/fopr.h"
#include "Measurements/FermionicM/source_types.hpp"
#include "Measurements/FermionicM/qprop_MultiShift.hpp"
#include "Dirac_ops/dirac_wilson.hpp"
#include "Dirac_ops/dirac_clover.hpp"
#include "Dirac_ops/dirac_DomainWall.hpp"

#include <assert.h>
#include <vector>

using namespace std;

int Test_RationalApprox::run(){
  CCIO::cout << "Starting Rational Approximation test" << std::endl;

  // Test standard constructor
  RationalApprox_params PsParameters(10, 10, 1, 2, 40, 0.05, 1.0);
  /*
  PsParameters.numerator_deg   = 10;
  PsParameters.denominator_deg = 10;
  
  PsParameters.exponent_num = 1;
  PsParameters.exponent_den = 2;

  PsParameters.gmp_remez_precision = 40;
  PsParameters.lambda_low          = 0.05;
  PsParameters.lambda_high         = 1.0;
  */

  RationalApprox TestApprox(PsParameters);

  // Test XML constructor
  RationalApprox TestXMLApprox(RA_node);

  // Test output
  // Reconstruct and test against pow
  double x_test = 0.5;
  double exponent = (double)PsParameters.exponent_num/(double)PsParameters.exponent_den;
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
  vector<int> spos(4,0);
  Source_local<Format::Format_F> Source(spos,
                                        CommonPrms::instance()->Nvol());

 
  //Dirac* Kernel = new Dirac_Wilson(0.01, &(Gfield_.data));
  //Dirac* Kernel = new Dirac_Clover(0.01, 1.0, &(Gfield_.data));

  int N5d   = 6;
  double M0 = -1.6;
  double c  = 1.0;
  double b  = 1.0;
  double mq = 0.01;
  vector<double> omega(N5d,1.0);
  Dirac_optimalDomainWall* Kernel = new Dirac_optimalDomainWall(b,c,M0,mq,omega,&(Gfield_.data));
  Source_local<Format::Format_F> Source5d(spos,
					  CommonPrms::instance()->Nvol()*N5d); 

   // Definition of the Solver
  int    Niter= 2000;
  double stop_cond = 1.0e-24;
  MultiShiftSolver* Solver = 
    new MultiShiftSolver_CG(new Fopr_DdagD(Kernel),
                            stop_cond,
                            Niter);
  

  // Solver test
  xqs.resize(Res.size());
  for (int i = 0; i < Res.size(); ++i)
    xqs.push_back(Field(CommonPrms::instance()->Nvol()*N5d));
  double residual;
  int Nconv;
  Solver->solve(xqs, Source5d.mksrc(0,0), Poles, residual, Nconv);

  //xqs contains the solutions of (M - Poles[i])x = b

  // Reconstruct solution (M^dag M)^(1/2)
  Field solution, solution2; 
  Field temp;

  solution = Source5d.mksrc(0,0);
  solution *= TestXMLApprox.Const();

  for (int i = 0; i < Poles.size(); ++i){
    temp = xqs[i];
    temp *= Res[i];
    solution += temp;
  }

  // Apply again (M^dag M)^(1/2)
  Solver->solve(xqs, solution, Poles, residual, Nconv);  

  solution2 = solution;
  solution2 *= TestXMLApprox.Const();

  for (int i = 0; i < Poles.size(); ++i){
    temp = xqs[i];
    temp *= Res[i];
    solution2 += temp;
  }  
  ////////////////////////////////
  // Check answer
  // Compare with (M^dag M)
  Field reference_sol, diff_field;
  
  temp = Kernel->mult(Source5d.mksrc(0,0));
  reference_sol = Kernel->mult_dag(temp);
  
  diff_field = reference_sol;
  diff_field -= solution2;

  CCIO::cout << ":::::::: Check answer -- diff (norm) = "<< diff_field.norm() <<"\n";

}
