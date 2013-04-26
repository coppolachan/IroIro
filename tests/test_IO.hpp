//------------------------------------------------------------------------
/*!
 * @file test_IO.hpp
 *
 * @brief Test for input/output functions
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------
#ifndef TEST_IO_H_
#define TEST_IO_H_

#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"


class Test_IO: public TestGeneral{
private:
  XML::node IOnode;
  GaugeField& Gfield_;
public:
  Test_IO(XML::node node,
	  GaugeField& Gfield):IOnode(node),
			      Gfield_(Gfield){}

  int run();  
};


#endif
