//------------------------------------------------------------------------
/*!
 * @file test_sunRepres.hpp
 *
 * @brief Test for sunRep classes
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------
#ifndef TEST_SUNREP_H_
#define TEST_SUNREP_H_

#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"


class Test_sunRep: public TestGeneral{
private:
  XML::node sunRepnode;
  GaugeField& Gfield_;
public:
  Test_sunRep(XML::node node,
	      GaugeField& Gfield):sunRepnode(node),
				  Gfield_(Gfield){}

  int run();  
};


#endif
