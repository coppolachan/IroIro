/*! @file test_HadronSpectrum_Nf2p1.hpp
    @brief Declaration of test code for meson & baryon spectrum
*/
#ifndef TEST_HADRONSPECTRUM_NF2P1_INCLUDED
#define TEST_HADRONSPECTRUM_NF2P1_INCLUDED

#include "include/common_code.hpp"
#include "tests/tests.hpp"
#include "Measurements/FermionicM/baryonCorrelator.hpp"

class Test_HadronSpectrum_Nf2p1: public TestGeneral{
private:
  XML::node node_;
  GaugeField conf_;
  std::string output_;
  void output_meson( std::ofstream&,const std::vector<double>&,std::string);
  void output_baryon(std::ofstream&,const correl_t&,const correl_t&,
		     std::string);
public:
  Test_HadronSpectrum_Nf2p1(XML::node node,const GaugeField& conf,
			    const RandNum&,std::string file)
    :node_(node),conf_(conf),output_(file){}
  int run();
};

#endif


