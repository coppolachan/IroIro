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

#include "Action/action_gauge.h"


int Test_HMC::run(XML::node node){
  std::cout << "Starting HMCrun" << std::endl;
  
  //Using factories
  
  RNG_Env::RNG = RNG_Env::createRNGfactory(node);
  /*
  Integrators::Integr = 
    Integrators::createIntegratorFactory(node, Gfield_.Format);
  */
  
  std::vector<int> multip(3);
  multip[0]= 1;
  multip[1]= 2;
  multip[2]= 2;
  
  Field* CommonField = new Field(Gfield_.Format.size());
  //these are to be created by the action factory
  Dirac* OpNf2    = new Dirac_Wilson(1.0/6.0,CommonField);
  Dirac* OpRatio1 = new Dirac_Wilson(1.0,CommonField);
  Dirac* OpRatio2 = new Dirac_Wilson(1.0/6.0,CommonField);
  
  ActionLevel al_1, al_2, al_3;
  Action* Gauge = new ActionGaugeWilson(6.0, 
					Gfield_.Format, 
					CommonField);
  al_1.push_back(Gauge);
 
  Solver* SolvNf2 = new Solver_CG(10e-12,
				  1000,
				  new Fopr_DdagD(OpNf2));


  Action* Nf2Action = new Action_Nf2(CommonField,
				     OpNf2,
				     SolvNf2);
  
  al_2.push_back(Nf2Action);
  Solver* SolvRatio1 = new Solver_CG(10e-12,
				     1000,
				     new Fopr_DdagD(OpRatio1));
  Solver* SolvRatio2 = new Solver_CG(10e-12,
				     1000,
				     new Fopr_DdagD(OpRatio2));
  
  Action* RatioAction = new Action_Nf2_ratio(CommonField,
					     OpRatio1, 
					     OpRatio2,
					     SolvRatio1,
					     SolvRatio2);
  al_3.push_back(RatioAction);
  ActionSet ASet;
  ASet.push_back(al_3);
  ASet.push_back(al_2);
  ASet.push_back(al_1);
  
  MDexec* Integrator = new MDexec_leapfrog(8,
					   3,
					   0.1,
					   ASet,
					   multip,
					   Gfield_.Format,
					   CommonField);

				      
  HMCgeneral hmc_general(node, *Integrator);  


  //Initialization
  //HMCgeneral hmc_general(node);

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
