//--------------------------------------------------------------------
/*! @file heatbath.hpp
 * @brief Declaration of classes for HeatBath update
 */
//--------------------------------------------------------------------
#ifndef HEATBATH_INCLUDED
#define HEATBATH_INCLUDED

#include "Tools/randNum_Factory.h"
#include "include/pugi_interface.h"
#include "include/common_fields.hpp"
#include "HBActions.hpp"
#include <memory>

struct HeatBathParams {
  int Nsweeps;
  int ThermalizationSteps;
  int StartingConfig;
  int SaveInterval; //Setting to 0 does not save configurations
  std::string Filename_prefix;

  HeatBathParams(pugi::xml_node node) {
    // ---- Default values
    SaveInterval = 1;  
    ThermalizationSteps = 0;
    StartingConfig = 1;
    Filename_prefix = "Conf_";
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
  }  

};


class HeatBathGeneral{
private:
  std::auto_ptr<const RandNum> rand_;
  std::auto_ptr<HBAction> Action_;
  const HeatBathParams Params_;

public:
  HeatBathGeneral(pugi::xml_node node,HBAction& Act, const RandNum& rand_num)
    :rand_(&rand_num),
     Action_(&Act),
     Params_(HeatBathParams(node)){}

  HeatBathGeneral(pugi::xml_node node)
    :rand_(RNG_Env::RNG->getRandomNumGenerator()),
     Params_(HeatBathParams(node)){}

  ~HeatBathGeneral(){}

  void update(GaugeField& U)const;



};


#endif
