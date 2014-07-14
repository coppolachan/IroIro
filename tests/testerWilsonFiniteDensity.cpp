//------------------------------------------------------------------------
/*!
 * @file testerWilsonFiniteDensity.cpp 
 * @brief Main source code for testing the Dirac_Wilson_FiniteDensity class
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 Time-stamp: <2014-05-24 10:52:55 noaki>
 */
//------------------------------------------------------------------------
#include "tests.hpp"
#include "test_Wilson_FiniteDensity.hpp"

int main(int argc, char* argv[]){

  Test_Wilson_FiniteDensity WFD= TestEnv::StartUp<Test_Wilson_FiniteDensity>(argc, argv);
  WFD.run();

  /*
  CREATE_TEST(Test_Wilson_FiniteDensity);
  RUN_TEST;
  CLEAR_TEST; 
  */
  return 0;
}

