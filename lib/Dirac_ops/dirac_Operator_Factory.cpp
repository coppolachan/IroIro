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
      if (!strcmp(Dirac_name, "DiracWilson_Brillouin")) 
	return new DiracWilsonBrillouinFactory(node);
      if (!strcmp(Dirac_name, "DiracClover"))  
	return new DiracCloverFactory(node);
      if (!strcmp(Dirac_name, "DiracOptimalDomainWall5d"))  
	return new DiracDWF5dFactory(node);
      if (!strcmp(Dirac_name, "DiracOptimalDomainWall5dEvenOdd"))  
	return new DiracDWF5dEvenOddFactory(node);
      std::cerr<<"No Dirac Operator available with name ["
	       << Dirac_name << "]. Request by <" << node.name() << ">\n";
      abort();
    }else{
      std::cout<<"Mandatory node is missing in input file (Dirac Object)\n";
      abort();
    }
  }

  DiracDWF4dOperatorFactory* 
  createDiracDWF4dOperatorFactory(const XML::node node){
    if (node !=NULL) {
      const char* Dirac_name = node.attribute("name").value();
      if (!strcmp(Dirac_name, "DiracOptimalDomainWall4d"))  
	return new DiracDWF4DfullFactory(node);
      if (!strcmp(Dirac_name, "DiracOptimalDomainWall4d_eo"))  
	return new DiracDWF4DeoFactory(node);
      std::cerr<<"No Dirac Operator available with name ["
	       << Dirac_name << "]. Request by <" << node.name() << ">\n";
      abort();
    }else{
      std::cout<<"Mandatory node is missing in input file (Dirac Object)\n";
      abort();
    }
  }

  // this one forces to have a PauliVillars operator
  // useful in the case of the DomainWall5d action
  DiracDWF5dOperatorFactory* 
  createDiracDWF5dOperatorFactory(const XML::node node){
    if (node !=NULL) {
      const char* Dirac_name = node.attribute("name").value();

      if (!strcmp(Dirac_name, "")) {
	std::cerr << "No name provided for DWF Dirac Operator. Request by <"
		  << node.name() << ">\n";
	abort();
      }
      if (!strcmp(Dirac_name, "DiracOptimalDomainWall5d"))  
	return new DiracDWF5dFactory(node);
      if (!strcmp(Dirac_name, "DiracOptimalDomainWall5dEvenOdd"))  
	return new DiracDWF5dEvenOddFactory(node);
      
      std::cerr <<"No Dirac Operator available with name ["
		<< Dirac_name << "]. Request by <" << node.name() << ">\n";
      abort();
    }else{
      std::cout<<"Mandatory node is missing in input file (DWF Dirac Object)\n";
      abort();
    }
  }
  
  // another clasification of the operator where inversion is not considered
  // like in the case of eigenmode calculation.
  DiracOperatorFactory* createGeneralDiracOperatorFactory(const XML::node node){
    // temporal hack
    return createGeneralDiracWilsonLikeOperatorFactory(node);
  }

  DiracWilsonLikeOperatorFactory* 
  createGeneralDiracWilsonLikeOperatorFactory(const XML::node node){
    if (node !=NULL) {
      const char* DWL_name = node.attribute("name").value();

      if (!strcmp(DWL_name, "")) {
	std::cerr << "No name provided for DiracWilsonLike Operator. Request by <"
		  << node.name() << ">\n";
	abort();
      }
      if (!strcmp(DWL_name, "DiracWilson"))  
	return new DiracWilsonFactory(node);
      if (!strcmp(DWL_name, "DiracClover"))  
	return new DiracCloverFactory(node);
      if (!strcmp(DWL_name, "DiracOptimalDomainWall5d"))  
	return new DiracDWF5dFactory(node);
      if (!strcmp(DWL_name, "DiracOptimalDomainWall4d"))  
	return new DiracDWF4DfullFactory(node);
      std::cerr<<"No DiracWilsonLike Operator available with name ["
	       << DWL_name << "]. Request by <" << node.name() << ">\n";
      abort();
    }else{
      std::cout<<"Mandatory node is missing in input file (Dirac Object)\n";
      abort();
    }
  }

  DiracStaggeredEvenOddLikeOperatorFactory* 
  createDiracStaggeredEvenOddLikeOperatorFactory(const XML::node node){
    if (node !=NULL) {
      const char* Dirac_name = node.attribute("name").value();

      if (!strcmp(Dirac_name, "")) {
	std::cerr<<"No name provided for DiracStaggeredLike Operator. Request by <"
		 << node.name() << ">\n";
	abort();
      }
      if (!strcmp(Dirac_name, "DiracStaggeredEvenOdd"))  
	return new DiracStaggeredEvenOddFactory(node);
      abort();
    }else{
      std::cout<<"Mandatory node is missing in input file (DiracStaggered Object)\n";
      abort();
    }
  }

} //end of namespace
