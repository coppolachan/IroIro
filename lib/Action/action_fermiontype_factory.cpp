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
      if (!strcmp(Action_name, "FourFlavorStaggered")) 
	return new FourFlavorStaggeredActionFactory(node);
      if (!strcmp(Action_name, "FourFlavorStaggeredRatio")) 
	return new FourFlavorStaggeredRatioActionFactory(node);
      if (!strcmp(Action_name, "TwoFlavorsDomainWall_5D")) 
	return new TwoFlavorDomainWall5dActionFactory(node);

      if (!strcmp(Action_name, "NfFlavorsDomainWall_5D")) 
	return new NfFlavorDomainWall5dActionFactory(node);

      //BGQ specific improved routines
#ifdef IBM_BGQ_WILSON
      if (!strcmp(Action_name, "TwoFlavorsDomainWall_5D-EO_BGQ"))
	return new TwoFlavorDomainWall5dEO_BGQ_ActionFactory(node);
      if (!strcmp(Action_name, "TwoFlavorsRatioDomainWall_5D-EO_BGQ"))
	return new TwoFlavorRatioDomainWall5dEO_BGQ_ActionFactory(node);
      if (!strcmp(Action_name, "NfFlavorsDomainWall_5D-EO_BGQ")) 
	return new NfFlavorDomainWall5d_EO_BGQ_ActionFactory(node);
#endif

      std::cerr << "No Fermionic Action available with name ["
		<< Action_name << "]. Request by <" << node.name() << ">\n";
      abort();
      
    } else {
      std::cout << "Mandatory node is missing in input file (Action Object)\n";
      abort();

    }
  }
}
