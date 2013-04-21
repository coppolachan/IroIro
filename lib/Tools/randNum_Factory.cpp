/*!
  @file randNum_Factory.cpp 
*/
#include "randNum_Factory.h"
#include "Communicator/comm_io.hpp"
#include <string.h>

namespace RNG_Env {
  RandomNumberCreator* createRNGfactory(XML::node node){
    XML::descend(node,"RandomNumberGen");
    if(node !=NULL){
      const char* RNG_name = node.attribute("name").value();
    
      if(!strcmp(RNG_name,"Mersenne Twister"))
	return new RandNum_MT19937_Creator(node);

      if(!strcmp(RNG_name,"DC Mersenne Twister"))
	return new RandNum_DCMT_Creator(node);

      CCIO::cerr<< "No Random Number Generator available with name ["
		<< RNG_name << "]" << std::endl;
      abort();
    }else{
      CCIO::cout<< "Mandatory node <RandomNumberGen> is missing in input file"
		<< std::endl;
      abort();
    }
  }
}
