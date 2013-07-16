/*! @file test_MesonSpectrum.hpp
    @brief Declaration of test code for meson spectrum
*/
#ifndef TEST_MESONSPECTRUM_SPMOM_INCLUDED
#define TEST_MESONSPECTRUM_SPMOM_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests.hpp"

class Test_MesonSpectrum_SpMom {
  Measurements::Input input_;
public:
  Test_MesonSpectrum_SpMom(Measurements::Input input):input_(input){}
  int run();
};

#endif


