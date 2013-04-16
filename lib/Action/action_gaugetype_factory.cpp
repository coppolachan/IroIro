/*!
 * @file action_gaugetype_factory.cpp 
 *
 * @brief Definition of createGaugeActionFactory function
 *
 * Time-stamp: <2013-04-16 14:28:15 neo>
 */

#include <string.h>
#include "action_gaugetype_factory.hpp"
#include "include/errors.hpp"

namespace GaugeAction {
  GaugeActionFactory* createGaugeActionFactory(XML::node node){
    std::ostringstream NoActionErr;
    if (node !=NULL) {

      const char* Action_name = node.attribute("name").value();

      if (!strcmp(Action_name, "")) {
        NoActionErr << "No name provided for Gauge Action. Check your XML input file\n";
        Errors::XMLerr(NoActionErr);
      }
      
      /////////////////////////////////////////////////     
      if (!strcmp(Action_name, "Wilson")) { 
	return new WilsonGaugeActionFactory(node);
      } 
      if (!strcmp(Action_name, "WilsonAdjoint")) { 
	return new WilsonGaugeAdjointActionFactory(node);
      } 
      if (!strcmp(Action_name, "Rectangle")) { 
	return new RectGaugeActionFactory(node);
      } 
      if (!strcmp(Action_name, "Iwasaki")) { 
	return new IwasakiGaugeActionFactory(node);
      } 
      if (!strcmp(Action_name, "Symanzik")) { 
	return new SymanzikGaugeActionFactory(node);
      } 
      if (!strcmp(Action_name, "DBW2")) { 
	return new DBW2GaugeActionFactory(node);
      } 
      /////////////////////////////////////////////////

      // If no action is found with provided name 
      // execution reaches this point
      NoActionErr << "No Gauge Action available with name ["
		  << Action_name << "]. Request by <" << node.name() << "> node";
      Errors::XMLerr(NoActionErr);
    } else {
      std::ostringstream NodeErr;
      NodeErr << "Mandatory node is missing in input file (Action Object)\n"
	      << "Check correct spelling of Gauge Action name.\n"
	      << "Refer to documentation for available names.";
      Errors::XMLerr(NodeErr);
    } 
  }
}
