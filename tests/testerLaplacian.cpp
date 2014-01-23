//------------------------------------------------------------------------
/*!
 * @file testerLapH.cpp 
 * @brief Main source code for testing the Laplacian class
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------
#include "tests.hpp"
#include "test_Laplacian.hpp"

int main(int argc, char* argv[]){

  CREATE_TEST(Test_LapH);
  RUN_TEST;
  CLEAR_TEST;  

  return 0;
}

