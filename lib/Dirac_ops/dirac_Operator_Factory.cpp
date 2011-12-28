#include "dirac_Operator_Factory.hpp"

#include <string.h>
#include "Communicator/comm_io.hpp"

namespace DiracOperators {
  DiracWilsonLikeOperatorFactory* 
  createDiracWilsonLikeOperatorFactory(const XML::node node){
    if (node !=NULL) {
      
      const char* Dirac_name = node.attribute("name").value();

      if (!strcmp(Dirac_name, "")) {
	std::cerr << "No name provided for Dirac Operator. Request by <"
		  << node.name() << ">\n";
	abort();
      }
      
      if (!strcmp(Dirac_name, "DiracWilson"))  
	return new DiracWilsonFactory(node);
      if (!strcmp(Dirac_name, "DiracWilson_EvenOdd")) 
	return new DiracWilsonEvenOddFactory(node);
      if (!strcmp(Dirac_name, "DiracOptimalDomainWall4d"))  
	return new DiracDWF4dFactory(node);
      if (!strcmp(Dirac_name, "DiracOptimalDomainWall5d"))  
	return new DiracDWF5dFactory(node);
            
      std::cerr << "No Dirac Operator available with name ["
		<< Dirac_name << "]. Request by <" << node.name() << ">\n";
      abort();
      
    } else {
      std::cout << "Mandatory node is missing in input file (Dirac Object)\n";
      abort();
    }
  }
}
