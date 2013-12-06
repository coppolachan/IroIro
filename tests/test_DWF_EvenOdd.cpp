/*!
 * @file test_wilson_EvenOdd.cpp
 * @brief Tests for the propagators 
 */
#include <stdio.h>

#include "test_DWF_EvenOdd.hpp"
#include "Measurements/FermionicM/fermion_meas_factory.hpp"

#include "Dirac_ops/dirac_wilson_EvenOdd.hpp"
#include "Solver/solver_CG.hpp"
#include "Solver/solver_BiCGStab.hpp"
#include "Measurements/FermionicM/qprop_EvenOdd.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include "Tools/RandomNumGen/randNum_MT19937.h"
#include "Measurements/FermionicM/source_types.hpp"

using namespace std;
using namespace Format;

int Test_Wilson_EvenOdd::run(){

  ///// Generating source vector 

  // local source
  vector<int> spos(4,0); 
  //vector<int> spos(4,1); 
  Source_local<Format::Format_F> src(spos,CommonPrms::instance()->Nvol());

  /*
  // noise source
  unsigned long init[4]={0x123, 0x234, 0x345, 0x456};
  int length=4;
  RandNum_MT19937 rand(init, length);

  Source_wnoise<SiteIndex_eo,Format::Format_F> src(rand,CommonPrms::instance()->Nvol());
  */
  /*
  //wall source not working
  Source_wall<Format_F> src(0,CommonPrms::instance()->Nvol());
  */

  /************************************************************************************/
  Field* u = &(conf_.U);
  double M0 = -1.6;

  // creation of Dirac_Wilson operators 
  Dirac_Wilson Dw(   M0,u);
  Dirac_Wilson Dw_eo(M0,u,Dop::EOtag());
  Dirac_Wilson Dw_oe(M0,u,Dop::OEtag());

  int N5=4;
  double b=2.0;
  double c=0.0;
  double mq =0.05;

  std::vector<double> omega(N5,1.0);
  
  // creation of Domain-Wall Dirac operators 
  Dirac_DomainWall         Ddwf(   b,c,M0,mq,omega,&Dw,u);
  Dirac_DomainWall_EvenOdd Ddwf_eo(b,c,M0,mq,omega,&Dw_eo,&Dw_oe,u);

  SiteIndex_eo* ieo = SiteIndex_eo::instance();
  vector<int> esec= ieo->esec();
  vector<int> osec= ieo->osec();

  Format::Format_F fmt(CommonPrms::instance()->Nvol());
  valarray<size_t> ie = fmt.get_sub(esec);
  valarray<size_t> io = fmt.get_sub(osec);

  // source generation
  Field ff= src.mksrc(1,1);
  Field fe= src.mksrc(esec,1,1);
  Field fo= src.mksrc(osec,1,1);

  // test of spritting into even/odd
  /*
  for(int i=0;i<ie.size();++i)
    CCIO::cout<<" (fe["<<i<<"],ff["<<ie[i]<<"])=("<<fe[i]<<","<<ff[ie[i]]<<")"
	      <<" (fo["<<i<<"],ff["<<io[i]<<"])=("<<fo[i]<<","<<ff[io[i]]<<")\n";
  */

  // test of shift_field (up)
  /*
  Format::Format_F fmh(CommonPrms::instance()->Nvol()/2);

  ShiftField_up<Format::Format_F>      ffp(ff,&fmt,0); valarray<double> ffpv(ffp.getva());
  ShiftField_even_up<Format::Format_F> fep(fe,&fmh,0); valarray<double> fepv(fep.getva());
  ShiftField_odd_up<Format::Format_F>  fop(fo,&fmh,0); valarray<double> fopv(fop.getva());

  for(int i=0;i<ie.size();++i)
    CCIO::cout<<" (fep["<<i<<"],ffp["<<io[i]<<"])=("<<fepv[i]<<","<<ffpv[io[i]]<<")"
	      <<" (fop["<<i<<"],ffp["<<ie[i]<<"])=("<<fopv[i]<<","<<ffpv[ie[i]]<<")"
	      <<std::endl;  
  */

  // test of shift_field (dn)
  /*
  ShiftField_dn<Format::Format_F>      ffm(ff,&fmt,0); valarray<double> ffmv(ffm.getva());
  ShiftField_even_dn<Format::Format_F> fem(fe,&fmh,0); valarray<double> femv(fem.getva());
  ShiftField_odd_dn<Format::Format_F>  fom(fo,&fmh,0); valarray<double> fomv(fom.getva());

  for(int i=0;i<ie.size();++i)
    CCIO::cout<<" (fem["<<i<<"],ffm["<<io[i]<<"])=("<<femv[i]<<","<<ffmv[io[i]]<<")"
	      <<" (fom["<<i<<"],ffm["<<ie[i]<<"])=("<<fomv[i]<<","<<ffmv[ie[i]]<<")"
	      <<std::endl;  
  */

  // test of Dirac mult

  Field wf = Ddwf.mult(ff);
  Field we = Ddwf_eo.mult_ee(fe);  we+= Ddwf_eo.mult_eo(fo);  
  Field wo = Ddwf_eo.mult_oe(fe);  wo+= Ddwf_eo.mult_oo(fo);  

  for(int i=0;i<ie.size();++i)
    CCIO::cout<<" (we["<<i<<"],wf["<<ie[i]<<"])=("<<we[i]<<","<<wf[ie[i]]<<")"
	      <<" (wo["<<i<<"],wf["<<io[i]<<"])=("<<wo[i]<<","<<wf[io[i]]<<")\n";

  /*
  // test of Solver & Qprop
  int    Niter= 1000;
  double stop_cond = 1.0e-24;

  prop_t sq; // propagator
  Fopr_DdagD DdagD(&Deo);
  //Solver_BiCGStab solver(stop_cond,Niter,&DdagD);
  Solver_CG solver(stop_cond,Niter,&DdagD);  

  Qprop_EvenOdd qprop(&Deo,&solver);
  qprop.calc(sq,src);

  prop_t sq_full; // propagator
  Fopr_DdagD DdagDf(&D);
  Solver_CG solver_full(stop_cond,Niter,&DdagDf);  
  Qprop qprop_full(&D,&solver_full);
  qprop_full.calc(sq_full,src);

  CCIO::cout<<"quark propagator obtained"<<std::endl;
  
  // meson correlators
  GammaMatrices::Unit Gamma;//pion
  MesonCorrelator meson(Gamma,Gamma);

  vector<double> mcorr = meson.calculate<Format::Format_F>(sq,sq);  
  vector<double>::const_iterator it=mcorr.begin();
  int t=0;
  while(it!=mcorr.end()) CommunicatorItems::pprintf ("%d %.8e\n",t++, *it++);

  mcorr = meson.calculate<Format::Format_F>(sq_full,sq_full);  
  it=mcorr.begin();
  t=0;
  while(it!=mcorr.end()) CommunicatorItems::pprintf ("%d %.8e\n",t++, *it++);
  */

  return 0;
}

