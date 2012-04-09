//------------------------------------------------------------------------
/*!
 * @file test_HMC_noFactories.cpp
 *
 * @brief run() function for HMCgeneral class test without factories
 *
 * @author Jun-Ichi Noaki
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include <stdio.h>
#include <time.h>

#include "HMC/hmcGeneral.hpp"
#include "test_HMC.hpp"
#include "Tools/randNum_Factory.h"
#include "Action/action_gauge_wilson.hpp"
#include "Action/action_gauge_rect.hpp"
#include "Action/action_Nf.hpp"
#include "Action/action_Nf_ratio.hpp"
#include "HMC/mdExec_leapfrog.hpp"
#include "Dirac_ops/dirac_wilson.hpp"
#include "Solver/solver_CG.hpp"
#include "Solver/rationalSolver.hpp"

int Test_HMC::run(){
  CCIO::cout << "Starting HMCrun\n";
  
  //Using factories just for RNG
  RNG_Env::RNG = RNG_Env::createRNGfactory(HMC_node);
  
  std::vector<int> multip(2);
  multip[0]= 1;
  multip[1]= 2;
  
  Action_Nf_params RationalParams;
  RationalParams.n_pseudof_ = 2;
  RationalParams.n_flav_    = 2;
  RationalParams.degree_.resize(3);
  RationalParams.degree_[MetroStep] = 15;
  RationalParams.degree_[MDStep]    = 10;
  RationalParams.degree_[PFStep]    = 10;

  RationalParams.precision_.resize(3);
  RationalParams.precision_[MetroStep] = 60;
  RationalParams.precision_[MDStep]    = 40;
  RationalParams.precision_[PFStep]    = 40;

  RationalParams.b_low_.resize(3);
  RationalParams.b_low_[MetroStep] = 0.1;
  RationalParams.b_low_[MDStep]    = 0.1;
  RationalParams.b_low_[PFStep]    = 0.1;

  RationalParams.b_high_.resize(3);
  RationalParams.b_high_[MetroStep] = 2.0;
  RationalParams.b_high_[MDStep]    = 2.0;
  RationalParams.b_high_[PFStep]    = 2.0;



  GaugeField* CommonField = new GaugeField;

  Action* Gauge = new ActionGaugeWilson(5.0, CommonField);
  //Action* Gauge = new ActionGaugeRect(2.25, 3.648, -0.331, CommonField);

  DiracWilsonLike* Kernel   = new Dirac_Wilson(0.1,&(CommonField->data));

  
  ActionLevel al_1, al_2;
  al_1.push_back(Gauge);
   
  MultiShiftSolver* Solver = new MultiShiftSolver_CG(new Fopr_DdagD(Kernel),
						     1e-20,
						     1000);
  
  RationalSolver* SolverRational = new RationalSolver(Solver);
  
  Action* NfAction = new Action_Nf(CommonField,
				   Kernel,
				   SolverRational,
				   RationalParams);

  al_2.push_back(NfAction);
  
  ActionSet ASet;
  ASet.push_back(al_2);
  ASet.push_back(al_1);
  
  MDexec* Integrator = new MDexec_leapfrog(8,
					   10,
					   0.02,
					   ASet,
					   multip,
					   CommonField);
				      
  HMCgeneral hmc_general(HMC_node, *Integrator);  
  ////////////// HMC calculation /////////////////
  try{
    CCIO::cout<< "HMC starts\n";
    hmc_general.evolve(Gfield_);
  }catch(const char* error){
    CCIO::cerr << error << std::endl;
    return EXIT_FAILURE;
  }
  ////////////////////////////////////////////////  

  return 0;
}
