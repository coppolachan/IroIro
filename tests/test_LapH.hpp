/*!
 * @file test_LapH.hpp
 * @brief Test for LapH solver
 *
 * @author Guido Cossu
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 *
 * Time-stamp: <2014-01-31 16:20:16 neo>
 */

#ifndef TEST_LAPH_HPP_
#define TEST_LAPH_HPP_

#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"

class Test_LapH_Solver{
  const Measurements::Input input_;
public:
  Test_LapH_Solver(const Measurements::Input& input):input_(input){}
  int run();  
};



#endif
