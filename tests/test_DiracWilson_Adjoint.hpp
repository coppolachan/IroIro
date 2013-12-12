/*!
 * @file test_DiracWilson_Adjoint.hpp
 * @brief Test for Dirac_Wilson_Adjoint class
 * @author Jun-Ichi Noaki
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */

#ifndef TEST_DWA_H_
#define TEST_DWA_H_

#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"

class Test_DiracWilson_Adjoint{
  const Measurements::Input input_;
public:
  Test_DiracWilson_Adjoint(const Measurements::Input& input):input_(input){}
  int run();  
};


#endif
