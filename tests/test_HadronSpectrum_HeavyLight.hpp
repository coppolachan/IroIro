/*! @file test_HadronSpectrum_HeavyLight.hpp
    @brief Declaration of test code for meson & baryon spectrum
    using multi source smearing and multi mass
*/
#ifndef TEST_HADRONSPECTRUM_HEAVYLIGHT_INCLUDED
#define TEST_HADRONSPECTRUM_HEAVYLIGHT_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"

class Test_HadronSpectrum_HeavyLight{
  const Measurements::Input input_;
public:
  Test_HadronSpectrum_HeavyLight(const Measurements::Input& input)
    :input_(input){}
  int run();
};

#endif


