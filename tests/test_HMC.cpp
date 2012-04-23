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

#include "HMC/hmcGeneral.hpp"
#include "test_HMC.hpp"

#include "Smearing/smearingFactories.hpp"

int Test_HMC::run(){
  CCIO::cout << "Starting HMCrun" << std::endl;
  
  RNG_Env::RNG = RNG_Env::createRNGfactory(HMC_node);

  Integrators::Integr = Integrators::createIntegratorFactory(HMC_node);

  //Initialization
  HMCgeneral hmc_general(HMC_node);

  ////////////// HMC calculation /////////////////
  double elapsed_time;
  TIMING_START;
  try{
    CCIO::cout<< "HMC starts"<<std::endl;
    hmc_general.evolve(Gfield_);
  }catch(const char* error){
    CCIO::cerr << error << std::endl;
    return EXIT_FAILURE;
  }

  TIMING_END(elapsed_time);
  CCIO::cout << "Total elapsed time (s): "<< elapsed_time/1000.0 << "\n";

  CCIO::cout << "Saving configuration on disk in binary format\n";
  CCIO::SaveOnDisk< Format::Format_G >(Gfield_.data, "final_conf.bin");
  return 0;
}
