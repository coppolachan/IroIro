#include "solver_Factory.hpp"

#include <string.h>

namespace SolverOperators {
SolverOperatorCreator* createSolverOperatorFactory(const XML::node node)
  {
    if (node !=NULL) {
      
      const char* Solver_name = node.attribute("type").value();
      
      if (!strcmp(Solver_name, "Solver_CG")) { 
	return new SolverCGCreator(node);
      }
      if (!strcmp(Solver_name, "Solver_BiCGStab")) { 
	return new SolverBiCGStabCreator(node);
      }
      std::cerr << "No Solver available with type ["
                << Solver_name << "]\n" 
		<< "Request by " << node.parent().name() << std::endl;
      abort();
    
    }else {
      std::cout << "Requested node is missing in input file (Solver Object)\n"
                << "Request by " << node.parent().name() << std::endl;
      abort();
    }
    
    
  };
  
};
