/*!
 * @file test_wilson_EvenOdd.cpp
 * @brief Tests for the propagators 
 */
#include <stdio.h>

#include "test_wilson_EvenOdd.hpp"
#include "Measurements/FermionicM/fermion_meas_factory_abs.hpp"
#include "Dirac_ops/dirac_wilson_EvenOdd.hpp"
#include "Solver/solver_CG.hpp"
#include "Solver/solver_BiCGStab.hpp"
#include "Measurements/FermionicM/qprop.hpp"
#include "Measurements/FermionicM/qprop_EvenOdd.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include "Tools/randNum_MT19937.h"
#include "Measurements/FermionicM/source_types.hpp"

using namespace std;
using namespace Format;

int Test_Wilson_EvenOdd::run(){

  // Generating source vector -------------------------------------------
  // local source
  /*
  vector<int> spos(4,0); 
  Source_local<Format::Format_F> src(spos,CommonPrms::instance()->Nvol());
  */
  // noise source
  unsigned long init[4]={0x123, 0x234, 0x345, 0x456};
  int length=4;
  RandNum_MT19937 rand(init, length);
  Source_wnoise<SiteIndex,Format::Format_F> 
    src(rand,CommonPrms::instance()->Nvol());    


  //Dirac Kernel definition ---------------------------------------------
  Dirac_Wilson_EvenOdd Deo(1.0/6.0, &(conf_.data));
  Dirac_Wilson         D(  1.0/6.0, &(conf_.data));
  SiteIndex_eo* ieo = SiteIndex_eo::instance();
  vector<int> esec= ieo->esec();
  vector<int> osec= ieo->osec();

  Format::Format_F fmt(CommonPrms::instance()->Nvol());
  valarray<size_t> ie = fmt.get_sub(esec);
  valarray<size_t> io = fmt.get_sub(osec);

  // test of Dirac mult
  Field ff= src.mksrc(1,1);
  Field fe= src.mksrc(esec,1,1);
  Field fo= src.mksrc(osec,1,1);

  Field wf = D.mult(ff);
  Field we = Deo.mult_ee(fe); we+= Deo.mult_eo(fo);  
  Field wo = Deo.mult_oe(fe); wo+= Deo.mult_oo(fo);  

  for(int i=0;i<ie.size();++i)
    CCIO::cout<<" (we["<<i<<"],wf["<<ie[i]<<"])=("<<we[i]<<","<<wf[ie[i]]<<")"
	      <<" (wo["<<i<<"],wf["<<io[i]<<"])=("<<wo[i]<<","<<wf[io[i]]<<")"
	      <<std::endl;  

  //------------ e/o preconditioned solver ------------------
  int    Niter= 1000;
  double stop_cond = 1.0e-24;

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
  Niter= 1000;
  stop_cond = 1.0e-24;

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

