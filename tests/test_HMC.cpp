//------------------------------------------------------------------------
/*!
 * @file test_HMC.cpp
 *
 * @brief run() function for HMCgeneral class test
 *
 * @author Jun-Ichi Noaki
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include <stdio.h>
#include <time.h>

#include "HMC/hmcGeneral.h"
#include "test_HMC.hpp"


int Test_HMC::run(XML::node node){
  std::cout << "Starting HMCrun" << std::endl;
  
  RNG_Env::RNG = RNG_Env::createRNGfactory(node);
  Integrators::Integr = 
    Integrators::createIntegratorFactory(node, Gfield_.Format);



  //Initialization
  HMCgeneral hmc_general(node);

  ////////////// HMC calculation /////////////////
  clock_t start_t = clock();
  int nodeid = Communicator::instance()->nodeid();

  try{
    if(nodeid==0) std::cout<< "HMC starts"<<std::endl;
    hmc_general.evolve(Gfield_.U);
  }catch(const char* error){
    if(nodeid==0) std::cerr << error << std::endl;
    return EXIT_FAILURE;
  }

  clock_t end_t = clock();
  if(nodeid==0) 
    std::cout << (double)(end_t -start_t)/CLOCKS_PER_SEC << std::endl;
  

  

  return 0;
}
