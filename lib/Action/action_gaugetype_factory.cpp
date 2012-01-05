#include "action_gaugetype_factory.hpp"

#include <string.h>

namespace GaugeAction {
  GaugeActionFactory* createGaugeActionFactory(XML::node node)
  {
    if (node !=NULL) {
      const char* Action_name = node.attribute("name").value();
      if (!strcmp(Action_name, "")) {
        std::cerr << "No name provided for Gauge Action. Check your xml file\n";
        abort();
      }
      
      /////////////////////////////////////////////////     
      if (!strcmp(Action_name, "Wilson")) { 
	return new WilsonGaugeActionFactory(node);
      } 
      if (!strcmp(Action_name, "Rectangle")) { 
	return new RectGaugeActionFactory(node);
      } 
      /////////////////////////////////////////////////
      std::cerr << "No Gauge Action available with name ["
                << Action_name << "]. Request by <" << node.name() << ">\n";
      abort();
      
    } else {
      std::cout << "Mandatory node is missing in input file (Action Object)\n";
      abort();
    } 


  }
  
}
