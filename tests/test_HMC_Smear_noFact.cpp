//------------------------------------------------------------------------
/*!
 * @file test_HMC.cpp
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

#include "Smearing/SmartConf.hpp"

int Test_HMC::run(XML::node node){
  std::cout << "Starting HMCrun" << std::endl;
  
  //Using factories just for RNG
  
  RNG_Env::RNG = RNG_Env::createRNGfactory(node);
  
  std::vector<int> multip(2);
  multip[0]= 1;
  multip[1]= 1;
  
  Field* CommonField = new Field(Gfield_.Format.size());
  // --------------------------------- Smearing
  Smear_APE BaseAPE(0.1);
  Smear_Stout AnalyticSmear(BaseAPE);
  int Nsmear = 1; // Smearing levels
  const bool nosmear = false;
  const bool dosmear = true;

  SmartConf ThinField;//empty for thin links
  SmartConf FatField(Nsmear, AnalyticSmear);

  /////////////////////////////////////////////


  ActionLevel al_1, al_2;
  // Gauge action
  Action* Gauge = new ActionGaugeWilson(6.2, 
					Gfield_.Format, 
					FatField.select_conf(nosmear));
  al_1.push_back(Gauge);

  // Fermionic action
  //  DiracWilsonLike* OpNf2    = new Dirac_Clover(1.0/6.0,1.0,FatField.select_conf(dosmear));
  DiracWilsonLike* OpNf2    = new Dirac_Wilson(1.0/6.0,FatField.select_conf(dosmear));
 
  Solver* SolvNf2 = new Solver_CG(1e-14,
				  1000,
				  new Fopr_DdagD(OpNf2));
  
  Action* Nf2Action = new Action_Nf2(FatField.select_conf(dosmear),
				     OpNf2,
				     SolvNf2,
				     true,
				     &FatField);
  al_2.push_back(Nf2Action);

  ActionSet ASet;
  ASet.push_back(al_2);
  ASet.push_back(al_1);
  
  MDexec* Integrator = new MDexec_leapfrog(8,
					   5,
					   0.02,
					   ASet,
					   multip,
					   Gfield_.Format,
					   &FatField);

				      
  HMCgeneral hmc_general(node, *Integrator);  

  // Note:
  // The line *U_=U in mdExec_leapfrog.cpp (init)
  // must be substituted by a call to SmartConf::set_GaugeField()


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
