#include "action_fermiontype_factory.hpp"

#include <string.h>

namespace FermionAction {
  FermionActionFactory* createFermionActionFactory(XML::node node){
    if (node !=NULL) {
      
      const char* Action_name = node.attribute("name").value();
      
      if (!strcmp(Action_name, "TwoFlavors")) 
	return new TwoFlavorActionFactory(node);
      if (!strcmp(Action_name, "TwoFlavorsRatio")) 
	return new TwoFlavorRatioActionFactory(node);
      if (!strcmp(Action_name, "TwoFlavorsDomainWall")) 
	return new TwoFlavorDomainWallActionFactory(node);
      if (!strcmp(Action_name, "TwoFlavorsEvenOdd")) 
	return new TwoFlavorDomainWallActionFactory(node);

      std::cerr << "No Fermion Action available with name ["
		<< Action_name << "]" << std::endl;
      abort();
      
    } else {
      std::cout << "Mandatory node is missing in input file (Action Object)"
		<< std::endl;
      abort();

    }
  }
}
