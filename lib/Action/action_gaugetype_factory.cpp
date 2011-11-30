#include "action_gaugetype_factory.hpp"

#include <string.h>

namespace GaugeAction {
  GaugeActionFactory* createGaugeActionFactory(XML::node node)
  {
    //XML::descend(node, "Action");

    if (node !=NULL) {
      
      const char* Action_name = node.attribute("name").value();
      
      if (!strcmp(Action_name, "Wilson")) { 
	return new WilsonGaugeActionFactory(node);
      } 

    } 


  };
  
}
