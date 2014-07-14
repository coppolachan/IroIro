#include "xmlUtilities.hpp"
//#include "pugi_interface.h"
#include "Communicator/comm_io.hpp"
#include <string.h>

using namespace std;

namespace XML{
  void nullCheck(const XML::node& node,string objName){
    if(node==NULL){
      std::cout<< "Mandatory node is missing in input file ("
	       << objName
	       << "obj)\n";
      abort();
    }else{
      const char* nodeName = node.attribute("name").value();
      if(!strcmp(nodeName,"")){
	CCIO::cerr<< "No name provided for "
		  <<  objName
		  << " c. Request by <"
		  << node.name()
		  <<">\n";
        abort();
      }
    }
  }

  void stopMsg(const XML::node& node,string objName){
    const char* nodeName = node.attribute("name").value();

    CCIO::cerr<<"No "
	      << objName
	      << " available with name ["
	      << nodeName 
	      << "]. Request by <"
	      << node.name()
	      <<">\n";
    abort();
  }
  
}

