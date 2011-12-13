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

#include "include/common_code.hpp"
#include "tests/tests.hpp"


class Test_IO: public TestGeneral{
private:
  GaugeField& Gfield_;
public:
  Test_IO(GaugeField& Gfield):Gfield_(Gfield){}

  int run(XML::node);  
};


#endif
