/*!
 * @file test_wilson.cpp
 *
 * @brief Tests for the propagators 
 *
 */
#include "test_wilson.hpp"
#include "Measurements/FermionicM/fermion_meas_factory.hpp"


#include "Dirac_ops/dirac_wilson.h"
#include "Solver/solver_CG.h"
#include "Solver/solver_BiCGStab.h"
#include "Measurements/FermionicM/qprop.h"
#include "Measurements/FermionicM/source_types.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include "Tools/randNum_MT19937.h"
#include <stdio.h>

using namespace std;
using namespace Format;

int Test_Wilson::run(){
  prop_t sq;// propagator
  // Use factories to construct the propagator
  QuarkPropagator* QP;
  XML::descend(Wilson_node_, "QuarkPropagator");

  Staples Staple(conf_.Format);
  CCIO::cout << "Plaquette : " << Staple.plaquette(conf_.U) << std::endl;


  QuarkPropagatorFactory* QP_Factory = 
    QuarkPropagators::createQuarkPropagatorFactory(Wilson_node_);
  //////////////////////////////////////

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456};
  int length=4;
  RandNum_MT19937 rand(init, length);

  // source ????
  vector<int> spos(4,0); 
  //vector<int> spos(4,1); 

  Source_local<Format::Format_F> src(spos,CommonPrms::instance()->Nvol());
  //Source_wnoise<Format_F> src(rand,CommonPrms::instance()->Nvol());
  //Source_wall<Format_F> src(0,CommonPrms::instance()->Nvol());
  //wall source not working

  // Without factories -----------------------------------------------------
  // Dirac Kernel definition
  // Dirac* Kernel = new Dirac_Wilson(1.0/6.0, &(conf_.U));
  Dirac* Kernel = new Dirac_Clover(1.0/6.0, 1.0, &(conf_.U));
  Kernel->update_internal_state();

  // Solver definition
  int    Niter= 1000;
  double stop_cond = 1.0e-24;
  //Solver* SolverBiCGstab = new Solver_BiCGStab(stop_cond,
  //   					       Niter,
  //   					       new Fopr_DdagD(Kernel));
  
  Solver_CG* SolverCG = new Solver_CG(stop_cond,
				      Niter,
				      new Fopr_DdagD(Kernel));

  for(int s = 0; s < 4; ++s){
    for(int c = 0; c < 3; ++c){
      std::cout << "Force ["<<s<<","<<c<<"] "
		<<(Kernel->md_force(src.mksrc(s,c), Kernel->mult(src.mksrc(s,c))))[1]
		<< std::endl;
    }
  }
  // quark propagator
  // we force a type check on the Kernel (must be DdagD type).
 
  

  Qprop QuarkPropagator(Kernel,SolverCG);
  QuarkPropagator.calc(sq,src);
  //---------------------------------------------------------------------------
  
  //  QP = QP_Factory->getQuarkProp(conf_);
  //  QP->calc(sq,src);
  
  CCIO::cout<<"quark propagator obtained"<<std::endl;
  
  // meson correlators
  /*
  GammaMatrices::Unit Gamma;//pion
  MesonCorrelator meson(Gamma,Gamma);
  vector<double> mcorr = meson.calculate<Format::Format_F>(sq,sq);  
  vector<double>::const_iterator it=mcorr.begin();
  int t=0;
  while(it!=mcorr.end()) CommunicatorItems::pprintf ("%d %.8e\n",t++, *it++);
  */

  return 0;
}

