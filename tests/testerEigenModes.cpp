//------------------------------------------------------------------------
/*!
 * @file testerEigenModes.cpp 
 * @brief Main source code for testing the EigenModes classes
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "test_EigenModes_IRL.hpp"

int main(int argc, char* argv[]){

  CREATE_TEST(Test_EigenModes_IRL);
  RUN_TEST;
  CLEAR_TEST;

  return 0;
}

