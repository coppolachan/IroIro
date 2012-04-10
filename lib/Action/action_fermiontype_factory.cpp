#include "action_fermiontype_factory.hpp"

#include <string.h>

namespace FermionAction {
  FermionActionFactory* createFermionActionFactory(XML::node node){
    if (node !=NULL) {
      
      const char* Action_name = node.attribute("name").value();
      
      if (!strcmp(Action_name, "")) {
	std::cerr << "No name provided for Action. Check your xml file\n";
	abort();
      }
      if (!strcmp(Action_name, "TwoFlavors")) 
	return new TwoFlavorActionFactory(node);
      if (!strcmp(Action_name, "NfFlavors")) 
	return new NfFlavorsActionFactory(node);
      if (!strcmp(Action_name, "TwoFlavorsRatio")) 
	return new TwoFlavorRatioActionFactory(node);
      if (!strcmp(Action_name, "NfFlavorsRatio")) 
	return new NfFlavorRatioActionFactory(node);

      if (!strcmp(Action_name, "TwoFlavorsDomainWall_5D")) 
	return new TwoFlavorDomainWall5dActionFactory(node);
      if (!strcmp(Action_name, "NfFlavorsDomainWall_5D")) 
	return new NfFlavorDomainWall5dActionFactory(node);

      std::cerr << "No Fermionic Action available with name ["
		<< Action_name << "]. Request by <" << node.name() << ">\n";
      abort();
      
    } else {
      std::cout << "Mandatory node is missing in input file (Action Object)\n";
      abort();

    }
  }
}
