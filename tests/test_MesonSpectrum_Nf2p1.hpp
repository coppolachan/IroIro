/*! @file test_MesonSpectrum_Nf2p1.hpp
    @brief Declaration of test code for meson spectrum
*/
#ifndef TEST_MESONSPECTRUM_NF2P1_INCLUDED
#define TEST_MESONSPECTRUM_NF2P1_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"

class Test_MesonSpectrum_Nf2p1{
  const Measurements::Input input_;
public:
  Test_MesonSpectrum_Nf2p1(const Measurements::Input& input):input_(input){}
  int run();
};

#endif


