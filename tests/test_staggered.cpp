/*!
 * @file test_staggered.cpp
 * @brief Tests for the propagators 
 */
#include "test_staggered.hpp"
#include "Measurements/FermionicM/fermion_meas_factory_abs.hpp"

#include "Dirac_ops/dirac_staggered.hpp"
#include "Solver/solver_CG.hpp"
#include "Solver/solver_BiCGStab.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include "Measurements/FermionicM/qprop.hpp"
#include "Measurements/FermionicM/source_types.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include "Tools/randNum_MT19937.h"
#include <stdio.h>
#include <time.h>
using namespace std;
using namespace Format;

int Test_staggered::run(){
  prop_t sq;// propagator
  // Use factories to construct the propagator
  QuarkPropagator* QP;
  XML::descend(stagg_node_, "QuarkPropagator");

  Staples Staple;
  CCIO::cout << "Plaquette : " << Staple.plaquette(conf_) << std::endl;

   QuarkPropagatorFactory* QP_Factory = 
    QuarkPropagators::createQuarkPropagatorFactory(stagg_node_);
  //////////////////////////////////////

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456};
  int length=4;
  RandNum_MT19937 rand(init, length);

  // source 
  vector<int> spos(4,0); 

  Source_local<Format::Format_S> src(spos,CommonPrms::instance()->Nvol());
  //Source_wnoise<Format_S> src(rand,CommonPrms::instance()->Nvol());
  //Source_wall<Format_S> src(0,CommonPrms::instance()->Nvol());
  //wall source not working

  // Without factories -----------------------------------------------------
  // Dirac Kernel definition
  Dirac_staggered Kernel(1.0/6.0, &(conf_.data));

  // Solver definition
  int    Niter= 1000;
  double stop_cond = 1.0e-24;
  Solver_CG SolverCG(stop_cond,Niter,new Fopr_DdagD(Kernel));

  // quark propagator
  // we force a type check on the Kernel (must be DdagD type).
  Qprop QuarkPropagator(&Kernel,&SolverCG);
  QuarkPropagator.calc(sq,src);
  //---------------------------------------------------------------------------
  
  //  QP = QP_Factory->getQuarkProp(conf_);
  //  QP->calc(sq,src);
  
  CCIO::cout<<"quark propagator obtained"<<std::endl;
  
  // meson correlators
  MesonCorrelator meson(Pion);
  vector<double> mcorr = meson.calculate<Format::Format_S>(sq,sq);  
  vector<double>::const_iterator it=mcorr.begin();
  int t=0;
  while(it!=mcorr.end()) CCIO::cout << t++ << "  "<< *it++ << "\n";

  return 0;
}

