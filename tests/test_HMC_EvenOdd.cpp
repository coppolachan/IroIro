//------------------------------------------------------------------------
/*!
 * @file test_HMC_EvenOdd.cpp
 *
 * @brief run() function for HMC EvenOdd class test
 *
 * @author Jun-Ichi Noaki
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include <stdio.h>
#include <time.h>
#include <vector>

#include "HMC/hmcGeneral.hpp"
#include "test_HMC_EvenOdd.hpp"
#include "Action/action_gauge_wilson.hpp"
#include "Action/action_Nf2.hpp"
#include "Action/action_Nf2_ratio.hpp"
#include "Communicator/comm_io.hpp"

using namespace std;

int Test_HMC_EvenOdd::run(){
  CCIO::cout << "Starting HMCrun" << std::endl;
 
  RNG_Env::RNG = RNG_Env::createRNGfactory(HMC_node);
  
  std::vector<int> multip(3);
  multip[0]= 1;
  multip[1]= 2;
  multip[2]= 2;
  
  Field* CommonField = new Field(Gfield_.Format.size());

  double mq1 = 0.10;
  double mq2 = 0.20;

  DiracWilsonLike* D1 = new Dirac_Wilson_EvenOdd(mq1,CommonField);
  DiracWilsonLike* D2 = new Dirac_Wilson_EvenOdd(mq2,CommonField);
   
  // gauge term
  ActionLevel al_1, al_2, al_3;

  Action* Gauge 
    = new ActionGaugeWilson(6.0,Gfield_.Format,CommonField);
  al_1.push_back(Gauge);

  // pf1 term
  Solver* SolvNf2 = new Solver_CG(10e-12,1000, new Fopr_DdagD(D2));

  Action* Nf2Action = new Action_Nf2(CommonField,D2,SolvNf2);
  al_2.push_back(Nf2Action);

  // pf2 term
  Solver* SolvR1 = new Solver_CG(10e-12,1000,new Fopr_DdagD(D1));
  Solver* SolvR2 = new Solver_CG(10e-12,1000,new Fopr_DdagD(D2));

  Action* RatioAction = new Action_Nf2_ratio(CommonField,
					     D1,
					     D2,
					     SolvR1,
					     SolvR2);
  al_3.push_back(RatioAction);

  ActionSet ASet;
  ASet.push_back(al_3);
  ASet.push_back(al_2);
  ASet.push_back(al_1);
  
  MDexec* Integrator = new MDexec_leapfrog(8,
					   3,
					   0.01,
					   ASet,
					   multip,
					   Gfield_.Format,
					   CommonField);

  HMCgeneral hmc_general(HMC_node, *Integrator);  

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
  CCIO::cout << (double)(end_t -start_t)/CLOCKS_PER_SEC 
	     << std::endl;
  return 0;
}
