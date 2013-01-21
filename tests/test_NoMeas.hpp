/*! @file test_NoMeas.hpp
    @brief Declaration of test code which does NOTHING
*/
#ifndef TEST_NOMEAS_INCLUDED
#define TEST_NOMEAS_INCLUDED

#include "include/common_code.hpp"
#include "tests/tests.hpp"

class Test_NoMeas: public TestGeneral{
public:
  Test_NoMeas(XML::node,const GaugeField&,const RandNum&,std::string){}
  int run(){
    CCIO::cout<<"Nothing to measure.\n";
  }
};

#endif


