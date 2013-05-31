//------------------------------------------------------------------------
/*!
 * @file testerDWF.cpp 
 * @brief Main source code for testing the DomainWall classes
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 Time-stamp: <2013-05-03 22:10:06 noaki>
 */
//------------------------------------------------------------------------
#include "tests.hpp"
#include "test_DWF.hpp"

int main(int argc, char* argv[]){

  CREATE_TEST(Test_DWF);
  RUN_TEST;
  CLEAR_TEST;

  return 0;
}

