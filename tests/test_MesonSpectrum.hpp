/*! @file test_MesonSpectrum.hpp
    @brief Declaration of test code for meson spectrum
*/
#ifndef TEST_MESONSPECTRUM_INCLUDED
#define TEST_MESONSPECTRUM_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"

class Test_MesonSpectrum {
  Measurements::Input input_;
public:
  Test_MesonSpectrum(const Measurements::Input& input):input_(input){}
  int run();
};

#endif


