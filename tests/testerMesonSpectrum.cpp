/*! @file testerMesonSpectrum.cpp 
 * @brief Main source code for testing the MesonSpectrum classes
 */
#include "tests.hpp"
#include "test_MesonSpectrum.hpp"

int main(int argc, char* argv[]){
  
  CREATE_RUN_TEST(Test_MesonSpectrum);
  CLEAR_TEST;
  
  return 0;
}

/*
int main(int argc, char* argv[]){
  int ret = TestEnv::StartRun<Test_MesonSpectrum>(argc,argv);
  return ret;
}
*/
