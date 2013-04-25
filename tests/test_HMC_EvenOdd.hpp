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

#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"


class Test_HMC_EvenOdd: TestGeneral{
private:
  XML::node HMC_node;
  GaugeField& Gfield_;
public:
  Test_HMC_EvenOdd(XML::node node,
		   GaugeField& Gfield):HMC_node(node),
				       Gfield_(Gfield){}

  int run();  
};


#endif
