//------------------------------------------------------------------------
/*!
 * @file test_HeatBath.hpp
 *
 * @brief Test for HeatBath update classes
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------
#ifndef TEST_HB_H_
#define TEST_HB_H_

#include "include/iroiro_code.hpp"

class Test_HB{
private:
  const char* test_name;
  XML::node HB_node;
  GaugeField Gfield_;
public:
  Test_HB(XML::node node,GaugeField& Gfield):HB_node(node),
					     Gfield_(Gfield){
    test_name = "HB";
    XML::descend(HB_node, test_name);    
  }

  int run();  
};


#endif
