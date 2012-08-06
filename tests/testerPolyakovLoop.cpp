/*! @file testerPolyakovLoop.cpp 
 * @brief Main source code for testing the PolyakovLoop class
 */
#include "tests.hpp"
#include "test_PolyakovLoop.hpp"

int main(int argc, char* argv[]){
  
  CREATE_RUN_TEST(Test_PolyakovLoop);
  CLEAR_TEST;
  
  return 0;
}

/*
int main(int argc, char* argv[]){
  int ret = TestEnv::StartRun<Test_PolyakovLoop>(argc,argv);
  return ret;
}
*/
