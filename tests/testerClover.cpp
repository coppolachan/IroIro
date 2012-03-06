//------------------------------------------------------------------------
/*!
  @file testerClover.cpp 
  
  @brief Main source code for testing the Clover classes
  
  @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>

*/
//------------------------------------------------------------------------
#include "test_clover.hpp"

int main(int argc, char* argv[]){

  CREATE_TEST(Test_Clover);
  RUN_TEST;
  CLEAR_TEST;  
  
  return 0;
}

