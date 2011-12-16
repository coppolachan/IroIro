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
  CCIO::cout << "Starting HMCrun" << std::endl;
  
  RNG_Env::RNG = RNG_Env::createRNGfactory(node);
  Integrators::Integr = 
    Integrators::createIntegratorFactory(node, Gfield_.Format);

  //Initialization
  HMCgeneral hmc_general(node);

  ////////////// HMC calculation /////////////////
  clock_t start_t = clock();
  try{
    CCIO::cout<< "HMC starts"<<std::endl;
    hmc_general.evolve(Gfield_.U);
  }catch(const char* error){
    CCIO::cerr << error << std::endl;
    return EXIT_FAILURE;
  }

  clock_t end_t = clock();
  CCIO::cout << (double)(end_t -start_t)/CLOCKS_PER_SEC << std::endl;

  //temporary
  CCIO::cout << "Saving configuration on disk in binary format\n";
  CCIO::SaveOnDisk< Format::Format_G >(Gfield_.U, "final_conf.bin");
  std::string seedfile("seed_file");
  //RNG_Env::RNG->getRandomNumGenerator()->saveSeed(seedfile);


  return 0;
}
