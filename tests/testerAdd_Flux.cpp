//------------------------------------------------------------------------
/*!
 * @file testerAdd_Flux.cpp 
 * @brief Main source code for Adding a U(1) flux to the Gauge configuration
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "test_Add_Flux.hpp"


int main(int argc, char* argv[]){
 
  CREATE_TEST(Test_Add_Flux);
  RUN_TEST;
  CLEAR_TEST;
  
  return 0;
}

