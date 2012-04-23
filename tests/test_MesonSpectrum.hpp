/*! @file test_MesonSpectrum.hpp
    @brief Declaration of test code for meson spectrum
*/
#ifndef TEST_MESONSPECTRUM_INCLUDED
#define TEST_MESONSPECTRUM_INCLUDED

#include "include/common_code.hpp"
#include "tests/tests.hpp"

class Test_MesonSpectrum: public TestGeneral{
private:
  XML::node node_;
  GaugeField conf_;
  GaugeField smeared_u_;
  GaugeField fixed_u_;
public:
  Test_MesonSpectrum(XML::node node,GaugeField conf)
    :node_(node),conf_(conf),
     smeared_u_(conf.Nvol()),
     fixed_u_(conf.Nvol()){
    XML::descend(node_,"Measurement");
  }

  int run();
};

#endif


