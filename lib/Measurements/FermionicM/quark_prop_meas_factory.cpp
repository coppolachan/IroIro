/*!
 * @file quark_prop_meas_factory.cpp 
 * @brief QuarkPropagators measurements operators factories selector
 */
#include "quark_prop_meas_factory.hpp"
#include <string.h>

namespace QuarkPropagators {
  QuarkPropagatorFactory* createQuarkPropagatorFactory(XML::node node) {
    
    if (node !=NULL) {
      const char* QuarkProp_name = node.attribute("name").value();
      
      if (!strcmp(QuarkProp_name, "Qprop")) { 
        return new QPropFactory(node);
      }
      if (!strcmp(QuarkProp_name, "Qprop_EvenOdd")) { 
        return new QPropFactory_EvenOdd(node);
      }
      if (!strcmp(QuarkProp_name, "QpropDWF")) { 
        return new QPropDWFFactory(node);
      }
      std::cerr << "No Quark Propagator available with name ["
		<< QuarkProp_name << "]" << std::endl;
      abort();
    } else {
      std::cout << "Requested node is missing in input file "
		 << "(QuarkPropagator Object)\n" 
		 << "Request by " << node.parent().name() << std::endl;
      abort();
    }
  }  
}
