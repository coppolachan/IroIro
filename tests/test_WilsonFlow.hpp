/*! @file test_WilsonFlow.hpp
    @brief Declaration of test code for the wilson flow
*/
#ifndef TEST_WILSONFLOW_INCLUDED
#define TEST_WILSONFLOW_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests.hpp"

class Test_WilsonFlow{
  const Measurements::Input input_;
public:
  Test_WilsonFlow(const Measurements::Input& input):input_(input){}
  void monitor(const GaugeField&,std::vector<double> &) const;
  void topologyoutput(const std::string&,int,std::vector<double> &) const;
  int run();
};

#endif
