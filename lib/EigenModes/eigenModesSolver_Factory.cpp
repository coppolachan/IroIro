/*! @file eigenModesSolver_Factory.cpp                       
 *  @brief EigenModeSolver-factories selector 
 */
#include "eigenModesSolver_Factory.hpp"
#include <string.h>

namespace EigenSolver{
  EigenSolverFactory* createEigenSolverFactory(XML::node node){
    if(node != NULL){
      const char* eslv_name = node.attribute("name").value();
      
      if(!strcmp(eslv_name, "ImplicitRestartedLanczos"))
	return new EigenSolverFactory_IRL(node);

      CCIO::cerr<<"No EigenSolver available with name["
		<< eslv_name << "]" << std::endl;
      abort();
    }else{
      CCIO::cout << "Requested node is missing in input file "
                 << "(EigenSolver Object)\n"
                 << "Request by " << node.parent().name() << std::endl;
      abort();
    }
  }
}
