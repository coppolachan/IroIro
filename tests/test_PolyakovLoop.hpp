/*! @file test_Polyakov.hpp
 *  @brief Declaration of classes for testing the PolyakovLoop classes
 */
#ifndef TEST_POLYAKOVLOOP_INCLUDED
#define TEST_POLYAKOVLOOP_INCLUDED

#include "include/common_code.hpp"
#include "Communicator/comm_io.hpp"

class Test_PolyakovLoop{
 private:
  XML::node node_;
  GaugeField conf_;
  std::string output_;
 public:
  Test_PolyakovLoop(XML::node node,const GaugeField& conf,
		    const RandNum&, std::string file)
    :node_(node),conf_(conf),output_(file){
    CCIO::cout<<"Test_PolyakovLoop called"<<std::endl;
  }
  int run();
};

#endif
