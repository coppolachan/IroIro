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
#include <assert.h>
#include <vector>

using namespace std;

int Test_RationalApprox::run(){
  CCIO::cout << "Starting Rational Approximation test" << std::endl;

  // Test standard constructor
  RationalApprox_params PsParameters;

  PsParameters.numerator_deg   = 10;
  PsParameters.denominator_deg = 10;
  
  PsParameters.exponent_num = 1;
  PsParameters.exponent_den = 2;

  PsParameters.gmp_remez_precision = 40;
  PsParameters.lambda_low          = 0.05;
  PsParameters.lambda_high         = 1.0;

  RationalApprox TestApprox(PsParameters);



  // Test XML constructor


  // Test output
  // Reconstruct and test against pow
  double x_test = 0.5;
  double exponent = (double)PsParameters.exponent_num/(double)PsParameters.exponent_den;
  double reference = pow(x_test, exponent);
  
  CCIO::cout << "Reference = "<< reference << "\n";

  //Reconstruct rational expansion

  double result;
  vector<double> Res = TestApprox.Residuals();
  vector<double> Poles = TestApprox.Poles();
  assert(Res.size() == Poles.size());

  result += TestApprox.Const();
  
  for (int i = 0; i < Res.size(); ++i) {
    result += Res[i]/(x_test + Poles[i]);
  }

  CCIO::cout << "Result = "<< result << "\n";
  CCIO::cout << "Difference = "<< result-reference << "\n";
  

}
