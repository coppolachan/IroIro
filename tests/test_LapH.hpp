/*!
 * @file test_LapH.hpp
 * @brief Test for LapH solver
 *
 * @author Guido Cossu
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */

#ifndef TEST_LAPH_H_
#define TEST_LAPH_H_

#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"

class Test_LapH_Solver{
  const char* test_name;
  XML::node LapH_node;
  GaugeField Gfield_;
public:
  Test_LapH_Solver(XML::node node,GaugeField& Gfield):LapH_node(node),
                                              Gfield_(Gfield){
    test_name = "LapH";
    XML::descend(LapH_node, test_name);    
  }



  int run();  
};


#endif
