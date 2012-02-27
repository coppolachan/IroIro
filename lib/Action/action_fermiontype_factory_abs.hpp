/*!
 * @file action_fermiontype_factory_abs.hpp 
 *
 * @brief Declaration of FermionType action factories
 */
#ifndef ACTION_FERMION_FACT_ABS_
#define ACTION_FERMION_FACT_ABS_

#include "include/pugi_interface.h"
#include "action_Factory.hpp"

class FermionActionFactory : public ActionFactory {
  virtual Action* getFermionAction(GaugeField* const) = 0;
public:
  Action* getAction(GaugeField* const F) {
    return getFermionAction(F);
  }
};

////////////////////////////////////////////////////
namespace FermionAction {
  FermionActionFactory* createFermionActionFactory(XML::node);

}

#endif
