//------------------------------------------------------------------------
/*!
 * @file testerDWF.cpp 
 * @brief Main source code for testing the DomainWall classes
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "test_DomainWall.hpp"

int main(int argc, char* argv[]){

  CREATE_TEST(Test_optimalDomainWall);
  RUN_TEST;
  CLEAR_TEST;

  return 0;
}

