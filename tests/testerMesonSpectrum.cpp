/*! @file testerMesonSpectrum.cpp 
 * @brief Main source code for testing the MesonSpectrum classes
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */

#include "test_MesonSpectrum.hpp"

int main(int argc, char* argv[]){
 
  CREATE_RUN_TEST(Test_MesonSpectrum);
  CLEAR_TEST;
  
  return 0;
}

/*
int main(int argc, char* argv[]){
  TestEnv::StartRun<Test_MesonSpectrum>(argc,argv);
  return 0;
}
*/
