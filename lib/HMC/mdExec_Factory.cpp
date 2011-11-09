#include "mdExec_Factory.h"

#include <string.h>


namespace Integrators{
  MDIntegratorCreator* createIntegratorFactory(XML::node node,
					       Format::Format_G& F){
    XML::descend(node, "Integrator");
    if (node !=NULL) {
      
      const char* Int_name = node.attribute("name").value();
      
      if (!strcmp(Int_name, "leapfrog_multistep")) { 
	return  new MDIntegrator_LeapfrogCreator(node, F);
      }
      std::cerr << "No Integrator available with name ["
		<< Int_name << "]" << std::endl;
      abort();
      
    } else {
      std::cout << "Mandatory node <Integrator> is missing in input file"
		<< std::endl;
      abort();
    }
  };
};
