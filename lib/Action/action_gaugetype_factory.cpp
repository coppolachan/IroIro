#include "action_gaugetype_factory.h"

#include <string.h>

namespace GaugeAction {
  GaugeActionCreator* createGaugeActionFactory(XML::node node)
  {
    //XML::descend(node, "Action");

    if (node !=NULL) {
      
      const char* Action_name = node.attribute("name").value();
      
      if (!strcmp(Action_name, "Wilson")) { 
	return new WilsonGaugeActionCreator(node);
      } 

    } 


  };
  
}
