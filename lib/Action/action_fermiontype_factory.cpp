#include "action_fermiontype_factory.hpp"

#include <string.h>

namespace FermionAction {
  FermionActionFactory* createFermionActionFactory(XML::node node)
  {
    if (node !=NULL) {
      
      const char* Action_name = node.attribute("name").value();
      
      if (!strcmp(Action_name, "TwoFlavors")) { 
	return new TwoFlavorActionFactory(node);
      }
      if (!strcmp(Action_name, "TwoFlavorsRatio")) { 
	return new TwoFlavorRatioActionFactory(node);
      } 
      
    } 
    
    
  };
  
};
