//------------------------------------------------------------------------
/*!
 * @file testerHMC.cpp 
 * @brief Main source code for testing the %HMC class
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "test_HMC.hpp"

int main(int argc, char* argv[]){

  CREATE_TEST(Test_HMC);
  RUN_TEST;
  CLEAR_TEST;  

  return 0;
}

