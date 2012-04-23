/*! @file testerMesonSpectrum.cpp 
 * @brief Main source code for testing the MesonSpectrum classes
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */

#include "test_MesonSpectrum.hpp"

int main(int argc, char* argv[]){
 
  CREATE_TEST(Test_MesonSpectrum);
  RUN_TEST;
  CLEAR_TEST;
  
  return 0;
}

