/*!
 * @file action_fermiontype_factory_abs.hpp 
 *
 * @brief Declaration of the abstract fermion action factory
 *
 * Time-stamp: <2013-04-17 14:34:40 neo>
 */

#ifndef ACTION_FERMION_FACT_ABS_
#define ACTION_FERMION_FACT_ABS_

#include "action_Factory.hpp"

/*
 * @brief Abstract class for fermionic action creation
 */
class FermionActionFactory : public ActionFactory {
  virtual Action* getFermionAction(GaugeField* const, SmartConf* const) = 0;
public:
  Action* getAction(GaugeField* const F, SmartConf* const SC) {
    return getFermionAction(F, SC);
  }
};

#endif
