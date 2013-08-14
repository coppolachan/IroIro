/*!
  @file action_fermiontype_factory_creator.cpp

  @brief Definition of createFermionActionFactory function

  Time-stamp: <2013-08-14 09:52:29 cossu>
*/

#include "action_fermiontype_factory.hpp"
#include "include/errors.hpp"


namespace FermionAction {
  FermionActionFactory* createFermionActionFactory(XML::node node){
    std::ostringstream NoActionErr;
    if (node !=NULL) {
      
      const char* Action_name = node.attribute("name").value();
	
      if (!strcmp(Action_name, "")) {
	Errors::XMLerr("No name provided for Action. Check your XML input file");
      }

      //////////////////////////////////////////////////////////////////////
      if (!strcmp(Action_name, "TwoFlavors")) 
	return new TwoFlavorActionFactory(node);
	
      // RHMC
      if (!strcmp(Action_name, "NfFlavors")) 
	return new NfFlavorsActionFactory(node);
	
      if (!strcmp(Action_name, "TwoFlavorsRatio")) 
	return new TwoFlavorRatioActionFactory(node);
	
      // RHMC 
      if (!strcmp(Action_name, "NfFlavorsRatio")) 
	return new NfFlavorRatioActionFactory(node);
	
      // Staggered actions
      if (!strcmp(Action_name, "FourFlavorStaggered")) 
	return new FourFlavorStaggeredActionFactory(node);
      if (!strcmp(Action_name, "FourFlavorStaggeredRatio")) 
	return new FourFlavorStaggeredRatioActionFactory(node);
	
      // DomainWall specific actions
      if (!strcmp(Action_name, "TwoFlavorsDomainWall_5D")) 
	return new TwoFlavorDomainWall5dActionFactory(node);
      if (!strcmp(Action_name, "NfFlavorsDomainWall_5D")) 
	return new NfFlavorDomainWall5dActionFactory(node);
	
      //BGQ specific improved routines
#ifdef IBM_BGQ_WILSON
      if (!strcmp(Action_name, "TwoFlavorsDomainWall_5D-EO_BGQ"))
	return new TwoFlavorDomainWall5dEO_BGQ_ActionFactory(node);
#ifdef HAVE_LIBBFM
      if (!strcmp(Action_name, "TwoFlavorsDomainWall_5D-EO_BFM"))
	return new TwoFlavorDomainWall5dEO_BFM_ActionFactory(node);
#endif
      if (!strcmp(Action_name, "TwoFlavorsRatioDomainWall_5D-EO_BGQ"))
	return new TwoFlavorRatioDomainWall5dEO_BGQ_ActionFactory(node);
      if (!strcmp(Action_name, "NfFlavorsDomainWall_5D-EO_BGQ")) 
	return new NfFlavorDomainWall5d_EO_BGQ_ActionFactory(node);
#endif
      //////////////////////////////////////////////////////////////////////
	
      // If no action is found with provided name 
      // execution reaches this point
      NoActionErr << "No Fermionic Action available with name ["
		  << Action_name << "].\nRequest by <" << node.name() << "> node";
      Errors::XMLerr(NoActionErr);
    } else {
      std::ostringstream NodeErr;
      NodeErr << "Mandatory node is missing in input file (Action Object)\n"
	      << "Check correct spelling of Fermion Action name.\n"
	      << "Refer to documentation for available names.";
      Errors::XMLerr(NodeErr);
    }
  }
}
