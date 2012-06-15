/*! @file gfixFactory.cpp                                          
    @brief QuarkPropagators measurements operators factories selector         
 */
#include "gfixFactory.hpp"
#include <string.h>

namespace GaugeFix{

  GFixFactory* createGaugeFixingFactory(XML::node& node){

    if(node != NULL){
      const char* Gfix_name = node.attribute("name").value();
      
      if(!strcmp(Gfix_name,"NoFixing")){
	return new GFixFactory_Free;
      }
      if(!strcmp(Gfix_name,"Landau")){
	return new GFixFactory_Landau(node);
      }
      if(!strcmp(Gfix_name,"Coulomb")){
	return new GFixFactory_Coulomb(node);
      }
      CCIO::cerr << "No Gauge Fixing available with name ["
                << Gfix_name << "]" << std::endl;
    }else{
      CCIO::cout << "Requested node is missing in input file "
		 << "(Gauge Fixing Object)\n"
		 << "Request by " << node.parent().name() << std::endl;
      abort();
    }
  }
}
