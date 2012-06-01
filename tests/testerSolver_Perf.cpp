//------------------------------------------------------------------------
/*!
 * @file testerPerf.cpp 
 * @brief Main source code for testing the DomainWall classes performance
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "test_Solver_Perf.hpp"

int main(int argc, char* argv[]){

  CREATE_TEST(Test_Solver_Performance);
  RUN_TEST;
  CLEAR_TEST;

  return 0;
}

