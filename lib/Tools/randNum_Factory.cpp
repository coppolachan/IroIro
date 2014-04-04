/*!
  @file randNum_Factory.cpp 
*/
#include "randNum_Factory.hpp"
#include "Communicator/comm_io.hpp"
#include <string.h>

namespace RNG_Env {
  void RNGfactory::createRNGfactory(XML::node node){
    XML::descend(node,"RandomNumberGen");

    if(node !=NULL){
      const char* RNG_name = node.attribute("name").value();
    
      if(!strcmp(RNG_name,"Mersenne Twister")){
	RNG.save(new RandNum_MT19937_Creator(node));
	is_initialized = true;
	return;
      }

      #ifdef HAVE_LIBDCMT
      if(!strcmp(RNG_name,"DC Mersenne Twister")){
	RNG.save(new RandNum_DCMT_Creator(node));
	is_initialized = true;
	return;
      }
      #endif
      if(!strcmp(RNG_name,"None")) {
	is_initialized = true;
	RNG.save(new NoRNG());
	return;
      }

      CCIO::cerr<< "No Random Number Generator available with name ["
		<< RNG_name << "]" << std::endl;
      abort();
    }else{
      CCIO::cout<< "Mandatory node <RandomNumberGen> is missing in input file"
		<< std::endl;
      abort();
    }
  }

  void initialize(XML::node node){
    RandNumG::instance().createRNGfactory(node);
  }
  
  RandNum* RNGfactory::getRNG() {
    if (is_initialized) {
      return RNG.get()->getRandomNumGenerator();
    } else 
      {
	CCIO::cerr << "Fatal Error: Randon number generator was not initialized properly.\n";
	CCIO::cerr << "The call to RNG_Env::initialize(XML::node) is missing.\n";
	abort();
      }
  }

}
