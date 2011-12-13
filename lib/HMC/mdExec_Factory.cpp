#include "mdExec_Factory.hpp"
#include <string.h>

namespace Integrators{
  MDIntegratorFactory* createIntegratorFactory(XML::node node,
					       Format::Format_G& F){
    XML::descend(node, "Integrator");
    if(node !=NULL){
      const char* Int_name = node.attribute("name").value();
      
      if(!strcmp(Int_name, "leapfrog_multistep")) 
	return  new MDIntegrator_LeapfrogFactory(node, F);

      std::cerr << "No Integrator available with name ["
		<< Int_name << "]" << std::endl;
      abort();
      
    }else{
      std::cout << "Mandatory node <Integrator> is missing in input file"
		<< std::endl;
      abort();
    }
  }
}
