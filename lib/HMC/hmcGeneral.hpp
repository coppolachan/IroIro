//--------------------------------------------------------------------
/*! @file hmcGeneral.hpp
 *
 * @brief Declaration of classes for HMC update
 *
 */
//--------------------------------------------------------------------
#ifndef HMCEXEC_INCLUDED
#define HMCEXEC_INCLUDED

#include <memory>

#include "Tools/randNum_Factory.h"
#include "HMC/mdExec_Factory.hpp"
#include "include/pugi_interface.h"

class Field;
class MDexec;
class RandNum;

struct HMCGeneralParams {
  int Nsweeps;
  int ThermalizationSteps;
  int StartingConfig;
  int SaveInterval; //Setting to 0 does not save configurations
  std::string Filename_prefix;

  HMCGeneralParams(pugi::xml_node node) {
    // ---- Default values
    SaveInterval = 1;  
    ThermalizationSteps = 0;
    StartingConfig = 1;
    Filename_prefix = "Conf_";
    // -------------------
    XML::read(node, "Nsweeps", Nsweeps, MANDATORY);
    if(XML::read(node, "Thermalization", ThermalizationSteps))
      CCIO::cout << "Using default [Thermalization = "<< ThermalizationSteps << "]\n";
    if(XML::read(node, "StartingConfig", StartingConfig))
      CCIO::cout << "Using default [StartingConfig  = "<< StartingConfig << "]\n";
    if(XML::read(node, "SaveInterval", SaveInterval))
      CCIO::cout << "Using default [SaveInterval = "<< SaveInterval << "]\n";
    if(XML::read(node, "SavePrefix", Filename_prefix))
      CCIO::cout << "Using default [SavePrefix = "<< Filename_prefix << "]\n";
  }
};

class HMCgeneral{
private:
  std::auto_ptr<MDexec> md_;
  std::auto_ptr<const RandNum> rand_;
  const HMCGeneralParams Params;

  bool metropolis_test(const double Hdiff)const;
  double evolve_step(Field&)const;

public:
  HMCgeneral(pugi::xml_node node, 
	     MDexec& md_exec, 
	     const RandNum& rand_num)
    :md_(&md_exec),
     rand_(&rand_num),
     Params(HMCGeneralParams(node)){}
  
  HMCgeneral(pugi::xml_node node, 
	     MDexec& md_exec)
    :md_(&md_exec),
     rand_(RNG_Env::RNG->getRandomNumGenerator()),
     Params(HMCGeneralParams(node)){}
  
  HMCgeneral(pugi::xml_node node)
    :md_(Integrators::Integr->getMDIntegrator()),
     rand_(RNG_Env::RNG->getRandomNumGenerator()),
     Params(HMCGeneralParams(node)){}
  
  ~HMCgeneral(){}
  


  void evolve(Field& U)const;
};

#endif
