/*!@file inputConfig.hpp 
 * @brief definition of the class which contain input information on a specific config.
 */
#ifndef INPUTCONFIG_INCLUDED
#define INPUTCONFIG_INCLUDED
#include "common_fields.hpp"
#include "lib/EigenModes/eigenModes.hpp"

struct InputConfig{
  GaugeField* gconf;
  std::vector<EigenModes*> emodes;

  InputConfig():gconf(NULL){}
  InputConfig(GaugeField* gf):gconf(gf){}

  Field* getGconf()const{ return &(gconf->data);}  
  bool hasEmodes()const{ return emodes.size() > 0;}
};
#endif


