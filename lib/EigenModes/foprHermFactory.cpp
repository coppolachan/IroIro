/*! @file foprHermFactory.cpp
 *  @brief Fopr_Herm -selecter
 */
#include "include/foprHermFactory.hpp"

namespace FoprHerm{
  FoprHermFactory* createFoprHermFactory(XML::node node){
    if(node !=NULL){
      const char* fh_name = node.attribute("name").value();
      if (!strcmp(fh_name,"")) {
	std::cerr << "No name provided for Fopr_Herm. Request by <"
		  << node.name() << ">\n";
	abort();
      }
     if(!strcmp(fh_name,"Hx"))    return new FoprHermFactory_H(node);
     if(!strcmp(fh_name,"DdagD")) return new FoprHermFactory_DdagD(node); 
     if(!strcmp(fh_name,"DDdag")) return new FoprHermFactory_DDdag(node); 
     if(!strcmp(fh_name,"Chebyshev_DdagD")) 
       return new FoprHermFactory_Chebyshev_DdagD(node);
    }else{
      std::cout<<"Mandatory node is missing in input file (Fopr_Herm Object)\n";
      abort();
    }
  }
}

