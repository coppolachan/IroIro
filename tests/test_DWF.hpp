/*!
 * @file test_DWF.hpp
 * @brief Test for DWF classes
 *
 * @author Jun-Ichi Noaki
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */

#ifndef TEST_DWF_H_
#define TEST_DWF_H_

#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"

class Test_DWF{
  const Measurements::Input input_;
public:
  Test_DWF(const Measurements::Input& input):input_(input){}
  int run();  
};


#endif
