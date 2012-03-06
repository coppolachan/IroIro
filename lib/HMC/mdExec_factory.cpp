#include "mdExec_factory.hpp"
#include <string.h>

namespace Integrators{
  MDIntegratorFactory* createIntegratorFactory(XML::node node){
    XML::descend(node, "Integrator");
    if(node !=NULL){
      const char* Int_name = node.attribute("name").value();
      
      if(!strcmp(Int_name, "leapfrog_multistep")) 
	return  new MDIntegrator_LeapfrogFactory(node);

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
