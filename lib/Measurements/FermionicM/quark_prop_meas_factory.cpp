/*!
 * @file quark_prop_meas_factory.cpp 
 * @brief QuarkPropagators measurements operators factories selector
 */
#include "quark_prop_meas_factory.hpp"
#include <string.h>

namespace QuarkPropagators {
  QuarkPropagatorFactory* createQuarkPropagatorFactory(XML::node node) {
    
    if (node !=NULL) {
      const char* qprop_name = node.attribute("name").value();
      //CCIO::cout<<"qprop_name="<<qprop_name<<std::endl;
      
      if (!strcmp(qprop_name, "Qprop")) { 
        return new QPropFactory(node);
      }
      if (!strcmp(qprop_name, "Qprop_EvenOdd")) { 
        return new QPropFactory_EvenOdd(node);
      }
      if (!strcmp(qprop_name, "QpropDWF")) { 
        return new QPropDWFFactory(node);
      }
      CCIO::cerr << "No Quark Propagator available with name ["
		<< qprop_name << "]" << std::endl;
      abort();
    } else {
      CCIO::cout << "Requested node is missing in input file "
		 << "(QuarkPropagator Object)\n" 
		 << "Request by " << node.parent().name() << std::endl;
      abort();
    }
  }  
}
