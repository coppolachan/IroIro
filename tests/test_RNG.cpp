//------------------------------------------------------------------------
/*!
 * @file test_RNG.cpp 
 * @brief Main source code for Randon Number Generator tests
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------
#include "include/common_code.hpp"
#include "include/geometry.hpp"
#include "Tools/randNum_MT19937.h"
#include "Tools/RandomNumGen/dcmt_wrapper.hpp"
#include "Communicator/comm_io.hpp"

using namespace XML;

int main(){

  //Reading input file
  node top_node = getInputXML("test_RNG.xml");  

  //Initializing geometry using XML input
  Geometry geom(top_node);


  unsigned long int seeds[] = { 0x123, 0x234, 0x345, 0x456 };
  int length = 4;

  CCIO::cout << "Testing Mersenne Twister 19937\n";

  RandNum_MT19937* rand1 = new RandNum_MT19937(seeds, length);
  std::string seed_file = "seed_file";
  //RandNum_MT19937* rand1 = new RandNum_MT19937(seed_file);
  rand1->saveSeed(seed_file);

  RandNum_MT19937* rand2 = new RandNum_MT19937();
  rand2->loadSeed(seed_file);

  CCIO::cout << "Original RNG      From saved file     Difference" 
	     << std::endl;  
  for (int i = 0; i < 10; ++i){
    double r1, r2;
    r1 = rand1->get();
    r2 = rand2->get();
    CCIO::cout << r1 << " " 
	       << r2 << " " 
	       << r1-r2 << std::endl;
  }


  delete rand1;
  delete rand2;


  CCIO::cout << "Testing Dynamic Creation of Mersenne Twister\n";
  CCIO::cout << "Creation of the MT object\n";
  
  RandNum_DCMT* dcmt = new RandNum_DCMT(1234, 1232);

  for (int i = 0; i < 10; ++i){
    double r1;
    r1 = dcmt->get();
    CCIO::cout << r1 << std::endl;
  }
  

}
