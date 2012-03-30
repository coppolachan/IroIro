//------------------------------------------------------------------------
/*!
 * @file test_RationalApprox.hpp
 *
 * @brief Test for %RationalApprox class
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------
#ifndef TEST_RATIONALAPPROX_HPP_
#define TEST_RATIONALAPPROX_HPP_

#include "include/common_code.hpp"
#include "tests/tests.hpp"

class Test_RationalApprox: TestGeneral{
private:
  const char* test_name;
  XML::node RA_node;
  GaugeField Gfield_;
public:
  Test_RationalApprox(XML::node node,GaugeField& Gfield):RA_node(node),
							 Gfield_(Gfield){
    test_name = "RationalApprox";
    XML::descend(RA_node, test_name);    
  }
  
  int run();  
};


#endif
