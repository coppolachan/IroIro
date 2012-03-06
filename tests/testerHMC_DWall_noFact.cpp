//------------------------------------------------------------------------
/*!
 * @file testerHMC_DomainWall.cpp 
 * @brief Main source code for testing the %HMC Domain Wall classes without
 factories
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "test_HMC_DWall_noFact.hpp"
#include "include/commandline.hpp"

using namespace XML;
using namespace MapsEnv;

int main(int argc, char* argv[]){

  CREATE_TEST(Test_HMC_DomainWall);
  RUN_TEST;
  CLEAR_TEST;  

  return 0;
}
