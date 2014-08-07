/*! @file testerEigenModesSolver_Lap.cpp                                        
 *  @brief Main source code for testing the EigenModesSolver classes  
 * for Scalar operators such as Laplacian
 */
#include "tests.hpp"
#include "test_EigenModesSolver_Lap.hpp"

int main(int argc, char* argv[]){                                             
  //CREATE_RUN_TEST(Test_EigenModesSolver);
  //CLEAR_TEST;

  TestEnv::StartRun<Test_EigenModesSolver_Lap>(argc,argv);       
  return 0;                                                         
}  
