//------------------------------------------------------------------------
/*!
 * @file test_DWF.hpp
 * @brief Test for DWF classes
 *
 * @author Jun-Ichi Noaki
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------
#ifndef TEST_DWF_H_
#define TEST_DWF_H_

#include "include/iroiro_code.hpp"

class Test_DWF{
private:
  XML::node DWF_node;
  GaugeField Gfield_;
public:
  Test_DWF(XML::node node,GaugeField& Gfield):DWF_node(node),
					      Gfield_(Gfield){
    const char* test_name = "DWF";
    XML::descend(DWF_node,test_name);    
  }

  int run();  
};


#endif
