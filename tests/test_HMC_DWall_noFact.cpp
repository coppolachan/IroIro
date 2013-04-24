//------------------------------------------------------------------------
/*!
 * @file test_HMC_DWall_noFact.cpp
 *
 * @brief run() function for HMCgeneral class test (DWF case)
 *
 * @author Jun-Ichi Noaki
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include <stdio.h>
#include <time.h>
#include <vector>

#include "HMC/hmcGeneral.hpp"
#include "test_HMC_DWall_noFact.hpp"
#include "Action/action_gauge_wilson.hpp"
#include "Action/action_Nf2_ratio.hpp"
#include "Dirac_ops/dirac_DomainWall.hpp"
#include "Solver/solver_CG.hpp"
#include "HMC/mdExec_leapfrog.hpp"

using namespace std;

int Test_HMC_DomainWall::run(){
  CCIO::cout << "Starting HMCrun" << std::endl;
 
  RNG_Env::RNG = RNG_Env::createRNGfactory(HMC_DW_node);
  
  std::vector<int> multip(3);
  multip[0]= 1;
  multip[1]= 2;
  multip[2]= 2;
  
  GaugeField* CommonField = new GaugeField;

  int N5d = 6;
  double M0 = -1.8;
  double b = 2.0;
  double c = 0.0;
  double mq1 = 0.05;
  double mq2 = 0.10;
  vector<double> omega(N5d,1.0);

  const Field* u = &(CommonField->data); 
  Dirac_Wilson Dw(M0,u);
  Dirac_optimalDomainWall Ddwf1( b,c,M0,mq1,omega,&Dw,u);
  Dirac_optimalDomainWall DdwfPV(b,c,M0,1.0,omega,&Dw,u);
  Dirac_optimalDomainWall Ddwf2( b,c,M0,mq2,omega,&Dw,u);

  // gauge term
  ActionLevel al_1, al_2, al_3;
  Action* Gauge 
    = new ActionGaugeWilson(6.0,CommonField);
  al_1.push_back(Gauge);

  // pf1 term
  Fopr_DdagD DdagD2(&Ddwf2);
  Fopr_DdagD DdagD_PV(&DdwfPV);
  Solver_CG SolvNf2(10e-12,1000,&DdagD2);
  Solver_CG SolvNf2PV(10e-12,1000,&DdagD_PV);
  Action* Nf2Action 
    = new Action_Nf2_ratio(CommonField,&Ddwf2,&DdwfPV,&SolvNf2,&SolvNf2PV);
  al_2.push_back(Nf2Action);

  // pf2 term
  Fopr_DdagD DdagD1(&Ddwf1);
  Solver_CG SolvR1(10e-12,1000,&DdagD1);
  Solver_CG SolvR2(10e-12,1000,&DdagD2);
  Action* RatioAction 
    = new Action_Nf2_ratio(CommonField,&Ddwf1,&Ddwf2,&SolvR1,&SolvR2);
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
					   CommonField);
  

  HMCgeneral hmc_general(HMC_DW_node,*Integrator);  

  ////////////// HMC calculation /////////////////
  try{
    CCIO::cout<< "HMC starts"<<std::endl;
    hmc_general.evolve(Gfield_);
  }catch(const char* error){
    CCIO::cerr << error << std::endl;
    return EXIT_FAILURE;
  }
  /////////////////////////////////////////////////  
  delete CommonField;

  return 0;
}
