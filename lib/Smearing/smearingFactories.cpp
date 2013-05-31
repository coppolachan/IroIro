/*!
 * @file smearingFactories.cpp 
 * @brief Definition of Smearing Operators switcher
 */

#include "smearingFactories.hpp"

#include <string.h>

namespace Smearings {
  SmearingFactory* 
  createSmearingFactory(const XML::node node){
    if (node !=NULL) {
      const char* Smear_name = node.attribute("type").value();

      if (!strcmp(Smear_name, "")) {
	CCIO::cerr << "No name provided for Smearing Operator. Request by <"
		  << node.name() << ">\n";
	abort();
      }
      if(!strcmp(Smear_name,"APE"))   return new APESmearingFactory(node);
      if(!strcmp(Smear_name,"Stout")) return new StoutSmearingFactory(node);

      CCIO::cerr<<"No Smearing Operator available with name ["
	       << Smear_name << "]. Request by <" << node.name() << ">\n";
      abort();
    }else{
      CCIO::cout<<"Mandatory node is missing in input file (Smear Object)\n";
      abort();
    }
  }

}
