/*!
 * @file test_wilson_EvenOdd.cpp
 * @brief Tests for the propagators 
 */
#include <stdio.h>

#include "test_wilson_EvenOdd.hpp"
#include "Measurements/FermionicM/fermion_meas_factory.hpp"
#include "Dirac_ops/dirac_wilson_EvenOdd.h"
#include "Solver/solver_CG.h"
#include "Solver/solver_BiCGStab.h"
#include "Measurements/FermionicM/qprop_EvenOdd.h"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include "Tools/randNum_MT19937.h"
#include "Measurements/FermionicM/source_types.hpp"

using namespace std;
using namespace Format;

int Test_Wilson_EvenOdd::run(){

  // Generating source vector -------------------------------------------
  // local source
  vector<int> spos(4,0); 
  Source_local<Format::Format_F> src(spos,CommonPrms::instance()->Nvol());
  /*
  // noise source
  unsigned long init[4]={0x123, 0x234, 0x345, 0x456};
  int length=4;
  RandNum_MT19937 rand(init, length);
  Source_wnoise<SiteIndex_eo,Format::Format_F> 
  src(rand,CommonPrms::instance()->Nvol());    
  */

  //Dirac Kernel definition ---------------------------------------------
  Dirac_Wilson_EvenOdd Deo(1.0/6.0, &(conf_.U));
  Dirac_Wilson         D(  1.0/6.0, &(conf_.U));
  int    Niter= 1000;
  double stop_cond = 1.0e-24;
  
  //------------ e/o preconditioned solver ------------------

  CCIO::cout<<"e/o-index: calc of quark propagator starts "<<std::endl;
  prop_t sq; // propagator
  Fopr_DdagD DdagD(&Deo);
  //Solver_BiCGStab solver(stop_cond,Niter,&DdagD);
  Solver_CG solver(stop_cond,Niter,&DdagD);  
  Qprop_EvenOdd qprop(&Deo,&solver);
  qprop.calc(sq,src);
  CCIO::cout<<"e/o-index: quark propagator obtained "<<std::endl;

  // meson correlators
  MesonCorrelator meson(Pion);
  vector<double> mcorr = meson.calculate<Format::Format_F>(sq,sq);  
  vector<double>::const_iterator it=mcorr.begin();
  int t=0;
  while(it!=mcorr.end()) CCIO::cout<<t++<<" "<< *it++ << std::endl;

  //---------------- full-index solver ------------------

  CCIO::cout<<"full-index: calc of quark propagator starts "<<std::endl;
  prop_t sq_full; // propagator
  Fopr_DdagD DdagDf(&D);
  Solver_CG solver_full(stop_cond,Niter,&DdagDf);  
  Qprop qprop_full(&D,&solver_full);
  qprop_full.calc(sq_full,src);
  CCIO::cout<<"full-index: quark propagator obtained "<<std::endl;

  // meson correlators
  mcorr = meson.calculate<Format::Format_F>(sq_full,sq_full);  
  it=mcorr.begin();
  t=0;
  while(it!=mcorr.end()) CCIO::cout<<t++<<" "<< *it++ << std::endl;

  return 0;
}

