/*!
 * @file test_ReadEigenModes.hpp
 * @brief Reads eigenmodes
 *
 * @author Guido Cossu
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 *
 * Time-stamp: <2014-05-09 15:18:06 neo>
 */

#ifndef TEST_READEIG_HPP_
#define TEST_READEIG_HPP_

#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"

class Test_ReadEigenModes{
  const Measurements::Input input_;
public:
  Test_ReadEigenModes(const Measurements::Input& input):input_(input){}
  int run();  
};



#endif
