/*! @file test_HadronSpectrum_Nf2p1.hpp
    @brief Declaration of test code for meson & baryon spectrum
*/
#ifndef TEST_HADRONSPECTRUM_NF2P1_INCLUDED
#define TEST_HADRONSPECTRUM_NF2P1_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"
#include "Measurements/FermionicM/baryonCorrelator.hpp"

class Test_HadronSpectrum_Nf2p1{
  const Measurements::Input input_;
  void output_meson( std::ofstream&,const std::vector<double>&,std::string);
  void output_baryon(std::ofstream&,const correl_t&,const correl_t&,
		     std::string);
public:
  Test_HadronSpectrum_Nf2p1(const Measurements::Input& input):input_(input){}
  int run();
};

#endif


