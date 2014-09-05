/*! @file testerQnumSuscept.cpp 
 * @brief Main source code to obtain building blocks of quark number susceptibility
 */
#include "tests.hpp"
#include "test_QnumSuscept.hpp"

int main(int argc, char* argv[]){
  
  CREATE_RUN_TEST(Test_QnumSuscept);
  CLEAR_TEST;
  
  return 0;
}

/*
int main(int argc, char* argv[]){
  int ret = TestEnv::StartRun<Test_QnumSuscept>(argc,argv);
  return ret;
}
*/
