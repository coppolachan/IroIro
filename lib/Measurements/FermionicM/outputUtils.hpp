#ifndef OUTPUTUTILS_INCLUDED
#define OUTPUTUTILS_INCLUDED
#include "Communicator/comm_io.hpp"
#include "meson_correlator.hpp"
#include "baryonCorrelator.hpp"

namespace Hadrons{

  void output_meson(std::ofstream& writer,
		    const std::vector<double>& Correl,
		    std::string msg);
  
  void output_baryon(std::ofstream& writer,
		     const correl_t& Upper,const correl_t& Lower,
		     std::string msg);
  
  void mesonProp(const prop_t& sq_ud,const prop_t& sq_s,std::ofstream& writer);
  void mesonPropGeneral(const prop_t& sq_ud,const prop_t& sq_s,std::ofstream& writer);
  void mesonExtraProp(const prop_t& sq_ud,const prop_t& sq_s,std::ofstream& writer);
  void baryonProp(const prop_t& sq_ud,const prop_t& sq_s,std::ofstream& writer);
}

#endif
