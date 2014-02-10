//--------------------------------------------------------------------
/*! @file hmcGeneral.hpp
 * @brief Declaration of classes for HMC update
 */
//--------------------------------------------------------------------
#ifndef HMCEXEC_INCLUDED
#define HMCEXEC_INCLUDED

#include "Tools/randNum_Factory.hpp"
#include "include/pugi_interface.h"
#include "HMC/mdExec_factory_abs.hpp"
#include "TrajectoryInfo.hpp"
#include <memory>

class MDexec;
class RandNum;

struct HMCGeneralParams {
  int Nsweeps; /* @brief Number of sweeps in this run */
  int TotalSweeps; /* @brief If provided, the total number of sweeps */
  int ThermalizationSteps;
  int StartingConfig;
  int SaveInterval; //Setting to 0 does not save configurations
  std::string Filename_prefix;

  HMCGeneralParams(pugi::xml_node node, TrajInfo* ExternalInfo) {
    CCIO::cout << "Creating HMCGeneralParams\n";
    TrajInfo Info;
    // ---- Default values
    if (ExternalInfo != NULL) {
      Info = *ExternalInfo;
    }

    SaveInterval = Info.SaveInterval;  
    ThermalizationSteps = Info.ThermalizationSteps;
    StartingConfig = Info.StartingConfig;
    Filename_prefix = Info.Filename_prefix;
 

    // -------------------
    XML::read(node,"Nsweeps", Nsweeps, MANDATORY);
    if(XML::read(node,"Thermalization", ThermalizationSteps))
      CCIO::cout<< "Using default [Thermalization = "
		<< ThermalizationSteps << "]\n";
    if(XML::read(node, "StartingConfig", StartingConfig))
      CCIO::cout<< "Using default [StartingConfig  = "
		<< StartingConfig << "]\n";
    if(XML::read(node, "SaveInterval", SaveInterval))
      CCIO::cout<< "Using default [SaveInterval = "
		<< SaveInterval << "]\n";
    if(XML::read(node, "SavePrefix", Filename_prefix))
      CCIO::cout<< "Using default [SavePrefix = "
		<< Filename_prefix << "]\n";
    TotalSweeps = StartingConfig+Nsweeps+1; //default
    if(XML::read(node,"TotalSweeps", TotalSweeps))
       CCIO::cout<< "Using default [TotalSweeps = "
		<< Nsweeps << "]\n";     

    Communicator::instance()->sync();
    CCIO::cout << std::flush;
    if (TotalSweeps <= StartingConfig){
      Nsweeps = 0;
      CCIO::cout<< "Run completed\n"; 
      exit(1);
    } else {
      int next_sweeps = TotalSweeps - StartingConfig+1;
      if (next_sweeps < Nsweeps) {
	Nsweeps = next_sweeps;
	CCIO::cout<< "Override [Nsweeps = "<< Nsweeps << "]\n"; 
      }    
    }

  }
};

class HMCgeneral{
private:
  std::auto_ptr<MDexec> md_;
  std::auto_ptr<const RandNum> rand_;
  const HMCGeneralParams Params;

  bool metropolis_test(const double Hdiff)const;
  double evolve_step(GaugeField&)const;
public:
  HMCgeneral(pugi::xml_node node,
	     MDexec& md_exec,
	     const RandNum& rand_num,
	     TrajInfo* Info = NULL)
    :md_(&md_exec),
     rand_(&rand_num),
     Params(HMCGeneralParams(node, Info)){}
  
  HMCgeneral(pugi::xml_node node,
	     MDexec& md_exec,
	     TrajInfo* Info = NULL)
    :md_(&md_exec),
     rand_(RNG_Env::RandNumG::instance().getRNG()),
     Params(HMCGeneralParams(node,Info)){}
  
  HMCgeneral(pugi::xml_node node, 
	     TrajInfo* Info = NULL)
    :md_(Integrators::Integr->getMDIntegrator()),
     rand_(RNG_Env::RandNumG::instance().getRNG()),
     Params(HMCGeneralParams(node, Info)){}
  
  ~HMCgeneral(){}

  void evolve(GaugeField& U)const;
};

#endif
