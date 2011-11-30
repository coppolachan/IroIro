//------------------------------------------------------------------------
/*!
 * @file test_HMC.hpp
 *
 * @brief Test for %HMC update classes
 *
 * @author Jun-Ichi Noaki
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------
#ifndef TEST_HMC_H_
#define TEST_HMC_H_

#include "include/common_code.hpp"

#include "tests/tests.hpp"


class Test_HMC: TestGeneral{
private:
  GaugeField& Gfield_;
public:
  Test_HMC(GaugeField& Gfield):Gfield_(Gfield){}

  int run(XML::node);  
};


#endif
