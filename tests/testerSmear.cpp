//------------------------------------------------------------------------
/*!
 * @file testerSmear.cpp 
 * @brief Main source code for testing the Smear classes
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "test_smear.hpp"

int main(int argc, char* argv[]){

  CREATE_TEST(Test_Smear);
  RUN_TEST;
  CLEAR_TEST;  
  
  return 0;
}

