/*! @file testerGWrelEigen.cpp 
 * @brief Main source code for testing the GWrelEigen class
 */
#include "tests.hpp"
#include "test_GWrelEigen.hpp"

int main(int argc, char* argv[]){
  
  CREATE_RUN_TEST(Test_GWrelEigen);
  CLEAR_TEST;
  
  return 0;
}

/*
int main(int argc, char* argv[]){
  int ret = TestEnv::StartRun<Test_GWrelEigen>(argc,argv);
  return ret;
}
*/
