/*! @file test_MesonSpectrum_Nf2p1.hpp
    @brief Declaration of test code for meson spectrum
*/
#ifndef TEST_MESONSPECTRUM_NF2P1_INCLUDED
#define TEST_MESONSPECTRUM_NF2P1_INCLUDED

#include "include/common_code.hpp"
#include "tests/tests.hpp"

class Test_MesonSpectrum_Nf2p1: public TestGeneral{
private:
  XML::node node_;
  GaugeField conf_;
  std::string output_;
public:
  Test_MesonSpectrum_Nf2p1(XML::node node,const GaugeField& conf,
			   const RandNum&,std::string file)
    :node_(node),conf_(conf),output_(file){}
  int run();
};

#endif


