/*!
 * @file test_CorrelateEigenModes.hpp
 * @brief Calculates the inner product of two eigenmodes, even from different files.
 *
 * @author Guido Cossu
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 *
 * Time-stamp: <2015-03-16 17:24:46 neo>
 */

#ifndef TEST_CORREIG_HPP_
#define TEST_CORREIG_HPP_

#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"

class Test_CorrelateEigenModes{
  const Measurements::Input input_;
public:
  Test_CorrelateEigenModes(const Measurements::Input& input):input_(input){}
  int run();  
};



#endif
