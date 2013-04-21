/*!
 * @file test_clover.cpp
 * @brief Tests for the propagators 
 */
#include "test_clover.hpp"
#include "Measurements/FermionicM/fermion_meas_factory_abs.hpp"

#include "Dirac_ops/dirac_clover.hpp"
#include "Solver/solver_CG.hpp"
#include "Solver/solver_BiCGStab.hpp"
#include "Measurements/FermionicM/qprop.hpp"
#include "Measurements/FermionicM/source_types.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include "Tools/RandomNumGen/randNum_MT19937.h"
#include <stdio.h>

using namespace std;
using namespace Format;

int Test_Clover::run(){
  prop_t sq;// propagator
  // Use factories to construct the propagator
  QuarkPropagator* QP;
  XML::descend(clover_node_, "QuarkPropagator");

  Staples Staple;
  CCIO::cout << "Plaquette : " << Staple.plaquette(conf_) << std::endl;


  QuarkPropagatorFactory* QP_Factory = 
    QuarkPropagators::createQuarkPropagatorFactory(clover_node_);
  //////////////////////////////////////

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456};
  int length=4;
  RandNum_MT19937 rand(init, length);

  // source 
  vector<int> spos(4,0); 

  Source_local<Format::Format_F> src(spos,CommonPrms::instance()->Nvol());
  //Source_wnoise<Format_F> src(rand,CommonPrms::instance()->Nvol());
  //Source_wall<Format_F> src(0,CommonPrms::instance()->Nvol());
  //wall source not working

  // Without factories -----------------------------------------------------
  // Dirac Kernel definition
  Dirac* Kernel = new Dirac_Clover(1.0/6.0, 1.0, &(conf_.data));
  //Kernel->update_internal_state();

  // Solver definition
  int    Niter= 1000;
  double stop_cond = 1.0e-24;
  //Solver* SolverBiCGstab = new Solver_BiCGStab(stop_cond,
  //   					       Niter,
  //   					       new Fopr_DdagD(Kernel));
  
  Solver_CG* SolverCG = new Solver_CG(stop_cond,
				      Niter,
				      new Fopr_DdagD(Kernel));

  // quark propagator
  // we force a type check on the Kernel (must be DdagD type).
  Qprop QuarkPropagator(Kernel,SolverCG);
  QuarkPropagator.calc(sq,src);
  //---------------------------------------------------------------------------
  
  //  QP = QP_Factory->getQuarkProp(conf_);
  //  QP->calc(sq,src);
  
  CCIO::cout<<"quark propagator obtained"<<std::endl;
  
  // meson correlators
  
  MesonCorrelator meson(Pion);
  vector<double> mcorr = meson.calculate<Format::Format_F>(sq,sq);  
  vector<double>::const_iterator it=mcorr.begin();
  int t=0;
  while(it!=mcorr.end()) CCIO::cout << t++ << "  "<< *it++ << "\n";
  

  return 0;
}

