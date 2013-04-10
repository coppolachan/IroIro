//------------------------------------------------------------------------
/*!
 * @file testerHeatBath.cpp 
 * @brief Main source code for testing the HeatBath class
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------
#include "tests.hpp"
#include "test_HeatBath.hpp"

int main(int argc, char* argv[]){

  CREATE_TEST(Test_HB);
  RUN_TEST;
  CLEAR_TEST;  

  return 0;
}

