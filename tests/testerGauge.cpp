//------------------------------------------------------------------------
/*!
 * @file testerGauge.cpp 
 * @brief Main source code for testing the Gauge configuration classes
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include "test_Gauge.hpp"


int main(int argc, char* argv[]){
 
  CREATE_TEST(Test_Gauge);
  RUN_TEST;
  CLEAR_TEST;
  
  return 0;
}

