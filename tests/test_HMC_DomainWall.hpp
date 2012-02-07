//------------------------------------------------------------------------
/*!
 * @file test_HMC_DomainWall.hpp
 *
 * @brief Test for %HMC with Domain Wall fermions update classes
 *
 * @author Jun-Ichi Noaki
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------
#ifndef TEST_HMC_DOMAINWALL
#define TEST_HMC_DOMAINWALL

#include "include/common_code.hpp"

#include "tests/tests.hpp"

class Test_HMC_DomainWall: TestGeneral{
private:
  XML::node HMC_DW_node;
  GaugeField& Gfield_;
public:
  Test_HMC_DomainWall(XML::node node,
		      GaugeField& Gfield):HMC_DW_node(node),
					  Gfield_(Gfield){}
  
  int run();  
};


#endif
