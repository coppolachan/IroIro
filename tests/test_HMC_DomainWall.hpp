//------------------------------------------------------------------------
/*!
 * @file test_HMC_DomainWall.hpp
 *
 * @brief Test for %HMC update classes
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
  GaugeField& Gfield_;
public:
  Test_HMC_DomainWall(GaugeField& Gfield)
    :Gfield_(Gfield){}

  int run(XML::node);  
};


#endif
