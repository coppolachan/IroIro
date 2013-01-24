/*! @file test_WilsonFlow.hpp
    @brief Declaration of test code for the wilson flow
*/
#ifndef TEST_WILSONFLOW_INCLUDED
#define TEST_WILSONFLOW_INCLUDED

#include "include/common_code.hpp"
#include "tests.hpp"

class Test_WilsonFlow: public TestGeneral{
private:
  XML::node node_;
  GaugeField conf_;
  std::string output_;
public:
  Test_WilsonFlow(XML::node node,const GaugeField& conf,
		  const RandNum&,std::string file)
    :node_(node),conf_(conf),output_(file){}
  void monitor(const GaugeField&, std::vector<double> &) const;
  void topologyoutput(std::string&, int, std::vector<double> &) const;
  int run();
};

#endif
