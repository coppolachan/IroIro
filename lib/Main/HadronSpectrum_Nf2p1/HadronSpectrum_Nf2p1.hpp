/*! @file HadronSpectrum_Nf2p1.hpp
    @brief Declaration of production code for meson & baryon spectrum
*/
#ifndef HADRONSPECTRUM_NF2P1_INCLUDED
#define HADRONSPECTRUM_NF2P1_INCLUDED

#include "include/iroiro_code.hpp"
#include "Measurements/FermionicM/baryonCorrelator.hpp"
#include "lib/Measurements/measGeneral.hpp"

class HadronSpectrum_Nf2p1{
  const Measurements::Input input_;
  void output_meson( std::ofstream&,const std::vector<double>&,std::string);
  void output_baryon(std::ofstream&,const correl_t&,const correl_t&,
		     std::string);
public:
  HadronSpectrum_Nf2p1(const Measurements::Input& input):input_(input){}
  int run();
};

#endif


