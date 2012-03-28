//------------------------------------------------------------------------
/*!
 * @file testerRationalApprox.cpp 
 * @brief Main source code for testing the %RationalApprox class
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "test_RationalApprox.hpp"

int main(int argc, char* argv[]){

  CREATE_TEST(Test_RationalApprox);
  RUN_TEST;
  CLEAR_TEST;  
  
  return 0;
}

