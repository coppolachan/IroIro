//------------------------------------------------------------------------
/*!
 * @file test_HB.cpp
 *
 * @brief run() function for HeatBathGeneral class test
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 *
 * Time-stamp: <>
 */
//------------------------------------------------------------------------


#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include "UpdatingAlg/HeatBath/heatbath.hpp"
#include "test_HeatBath.hpp"

#include "UpdatingAlg/HeatBath/HBActions_SUN.hpp"


int Test_HB::run(){
  CCIO::cout << "Starting Heat Bath run" << std::endl;
  
  RNG_Env::RNG = RNG_Env::createRNGfactory(HB_node);

  //Initialization
  //HeatBathGeneral hb_general(HB_node);
  HBAction_SUN SUNAction(2.3, &Gfield_);
  HeatBathGeneral hb_general(HB_node, SUNAction,*(RNG_Env::RNG->getRandomNumGenerator()));

  ////////////// HMC calculation /////////////////
  double elapsed_time;
  TIMING_START;
  try{
    CCIO::cout<< "HB starts"<<std::endl;
    hb_general.update(Gfield_);
  }catch(const char* error){
    CCIO::cerr << error << std::endl;
    return EXIT_FAILURE;
  }

  TIMING_END(elapsed_time);
  CCIO::cout << "Total elapsed time (s): "<< elapsed_time/1000.0 << "\n";

  
  sleep(2);// Hack to avoid concurrency problems
  /*
  CCIO::cout << "Saving configuration on disk in binary format\n";
  CCIO::SaveOnDisk< Format::Format_G >(Gfield_.data, "final_conf.bin");
  */

  return 0;
}
