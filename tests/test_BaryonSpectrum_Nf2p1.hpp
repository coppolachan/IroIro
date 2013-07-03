/*! @file test_BaryonSpectrum_Nf2p1.hpp
    @brief Declaration of test code for baryon spectrum
*/
#ifndef TEST_BARYONSPECTRUM_NF2P1_INCLUDED
#define TEST_BARYONSPECTRUM_NF2P1_INCLUDED

#include "Measurements/FermionicM/baryonCorrelator.hpp"
#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"

class Test_BaryonSpectrum_Nf2p1{
private:
  const Measurements::Input input_;
  void output(std::ofstream&,const correl_t&,const correl_t&,std::string msg);
public:
  Test_BaryonSpectrum_Nf2p1(const Measurements::Input& input):input_(input){}
  int run();
};

#endif


