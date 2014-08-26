//------------------------------------------------------------------------
/*!
 * @file testerSolver_HDCG.cpp 
 * @brief Main source code for testing the BFM classes, HDCG solver
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "test_Solver_HDCG.hpp"

int main(int argc, char* argv[]){

  CREATE_TEST(Test_Solver_HDCG);
  RUN_TEST;
  CLEAR_TEST;

  return 0;
}

