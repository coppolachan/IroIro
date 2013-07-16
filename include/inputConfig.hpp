/*!@file inputConfig.hpp 
 * @brief definition of the class which contain input information on a specific config.
 */
#ifndef INPUTCONFIG_INCLUDED
#define INPUTCONFIG_INCLUDED
#include "common_fields.hpp"
#include "lib/EigenModes/eigenModes.hpp"

struct InputConfig{
  GaugeField* const gconf;
  EigenModes* const emode;
  
  InputConfig(GaugeField* const conf,EigenModes* ems=NULL)
    :gconf(conf),emode(ems){}

  Field* getGconf(){ return &(gconf->data);}
  std::vector<Field>* getEvec(){ return &(emode->evecs_);}
  std::vector<double>* getEval(){ return &(emode->evals_);}
};
#endif


