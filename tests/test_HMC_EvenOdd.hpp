//------------------------------------------------------------------------
/*!
 * @file test_HMC_EvenOdd.hpp
 *
 * @brief Test for %HMC update classes
 *
 * @author Jun-Ichi Noaki
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------
#ifndef TEST_HMC_EVENODD_H_
#define TEST_HMC_EVENODD_H_

#include "include/common_code.hpp"
#include "tests/tests.hpp"


class Test_HMC_EvenOdd: TestGeneral{
private:
  GaugeField& Gfield_;
public:
  Test_HMC_EvenOdd(GaugeField& Gfield):Gfield_(Gfield){}

  int run(XML::node);  
};


#endif
