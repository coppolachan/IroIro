//------------------------------------------------------------------------
/*!
 * @file testerLapH.cpp 
 * @brief Main source code for testing the LapH classes
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 * 
 * Time-stamp: <2014-01-31 16:19:56 neo>
 */
//------------------------------------------------------------------------
#include "tests.hpp"
#include "test_LapH.hpp"

int main(int argc, char* argv[]){

  CREATE_RUN_TEST(Test_LapH_Solver);
  CLEAR_TEST;

  return 0;
}

