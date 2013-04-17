/*!
 * @file action_fermiontype_factory_creator.hpp 
 *
 * @brief Declaration of the FermionAction namespace for factories
 *
 * Time-stamp: <2013-04-17 14:33:49 neo>
 */

#ifndef ACTION_FERMION_FACT_CREATOR_
#define ACTION_FERMION_FACT_CREATOR_

#include "include/pugi_interface.h"
#include "action_fermiontype_factory_abs.hpp"

namespace FermionAction {
  FermionActionFactory* createFermionActionFactory(XML::node);
}


#endif
