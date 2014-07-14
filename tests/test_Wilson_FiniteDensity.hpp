/*!
 * @file test_Wilson_FiniteDensity.hpp
 * @brief Tests for the mult & mult_dag
 */
#ifndef TEST_WILSON_FINITEDENSITY_INCLUDED
#define TEST_WILSON_FINITEDENSITY_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests.hpp"
#include "Dirac_ops/dirac_wilson_FiniteDensity.hpp"
class Test_Wilson_FiniteDensity{
private:
  const char* test_name;
  XML::node node_;
  GaugeField conf_;
public:
  Test_Wilson_FiniteDensity(const XML::node node, GaugeField conf)
    :node_(node),
     conf_(conf){
    CCIO::cout<<"conf[0]="<<conf.data[0]<<"\n";
    CCIO::cout<<"conf.size()="<<conf.size()<<"\n";
    CCIO::cout<<"conf_[0]="<<conf_.data[0]<<"\n";
    CCIO::cout<<"conf_.size()="<<conf_.size()<<"\n";
    test_name="WilsonFiniteDensity";
    XML::descend(node_,test_name); 
  }

  void run();

};

#endif
