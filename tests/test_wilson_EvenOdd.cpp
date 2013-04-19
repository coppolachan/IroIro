/*!
 * @file test_wilson_EvenOdd.cpp
 * @brief Tests for the propagators 
 */
#include <stdio.h>
#include <algorithm>


#include "test_wilson_EvenOdd.hpp"
#include "Measurements/FermionicM/fermion_meas_factory_abs.hpp"
#include "Dirac_ops/dirac_wilson_EvenOdd.hpp"
#include "Solver/solver_CG.hpp"
#include "Solver/solver_BiCGStab.hpp"
#include "Measurements/FermionicM/qprop.hpp"
#include "Measurements/FermionicM/qprop_EvenOdd.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include "Tools/RandomNumGen/randNum_MT19937.h"
#include "Measurements/FermionicM/source_types.hpp"
#include "Main/Geometry/siteMap.hpp"

using namespace std;
using namespace Format;
using namespace SiteMap;

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
  SiteIndex_EvenOdd* ieo = SiteIndex_EvenOdd::instance();

  vector<int> esec= ieo->esec();
  vector<int> osec= ieo->osec();

  Format::Format_F fmt(CommonPrms::instance()->Nvol());
  valarray<size_t> ie = fmt.get_sub(esec);
  valarray<size_t> io = fmt.get_sub(osec);

  // generation of the source
  Field ff= src.mksrc(1,1);
  Field fe= src.mksrc(esec,1,1);
  Field fo= src.mksrc(osec,1,1);

  int Nx = CommonPrms::instance()->Nx();
  int Ny = CommonPrms::instance()->Ny();
  int Nz = CommonPrms::instance()->Nz();

  /*
  //////// test of manual shift (forward)////////////

  // full site
  vector<int> bdry_b,bdry_u,bulk_b,bulk_u;
  Field pff(ff.size());  
  // x=0 (sender)
  for(int n=0;n<shiftSite.slice_size(0,XDIR);++n)
    bdry_b.push_back(shiftSite.xslice(0,n,XDIR));
  // x=Nx-1 (receiver)
  for(int n=0;n<shiftSite.slice_size(Nx-1,XDIR);++n)
    bdry_u.push_back(shiftSite.xslice(Nx-1,n,XDIR));
  // transfer
  valarray<double> recv(fmt.Nin()*bdry_b.size());
  Communicator::instance()->transfer_fw(recv,ff[fmt.get_sub(bdry_b)],XDIR);
  pff.set(fmt.get_sub(bdry_u),recv);
  // bulk
  for(int x=0;x<Nx-1;++x)
    for(int n=0;n<shiftSite.slice_size(x,XDIR);++n)
      pff.set(fmt.islice(shiftSite.xslice(x,n,XDIR)),
	      ff[fmt.islice(shiftSite.xslice(x+1,n,XDIR))]);

  // even/odd site
  Format::Format_F fmh(CommonPrms::instance()->Nvol()/2);

  vector<int> ebdry_b,ebdry_u;   // o->e 
  vector<int> obdry_b,obdry_u;   // e->o 
  valarray<double> recv_e,recv_o;
  Field pfe(fe.size());  
  Field pfo(fo.size()); 
  //x=0(sender)
  for(int n=0;n<shiftSite_eo.slice_size(0,XDIR);++n)
    ebdry_b.push_back(shiftSite_eo.xslice(0,n,XDIR));
  for(int n=0;n<shiftSite_oe.slice_size(0,XDIR);++n)
    obdry_b.push_back(shiftSite_oe.xslice(0,n,XDIR));

  //x=Nx-1(receiver)
  for(int n=0;n<shiftSite_eo.slice_size(Nx-1,XDIR);++n)
    ebdry_u.push_back(shiftSite_eo.xslice(Nx-1,n,XDIR));
  for(int n=0;n<shiftSite_oe.slice_size(Nx-1,XDIR);++n)
    obdry_u.push_back(shiftSite_oe.xslice(Nx-1,n,XDIR));
  
  for(int i=0; i<ebdry_b.size();++i)
    CCIO::cout<<i
	      <<" ebdry_b= "<<ebdry_b[i]
	      <<" ebdry_u= "<<ebdry_u[i]
	      <<" obdry_b= "<<obdry_b[i]
	      <<" obdry_u= "<<obdry_u[i]<<std::endl;
  
  // transfer
  recv_e.resize(fmh.Nin()*ebdry_b.size());
  Communicator::instance()->transfer_fw(recv_e,fo[fmh.get_sub(obdry_b)],XDIR);
  pfe.set(fmh.get_sub(ebdry_u),recv_e);
  recv_o.resize(fmh.Nin()*obdry_b.size());
  Communicator::instance()->transfer_fw(recv_o,fe[fmh.get_sub(ebdry_b)],XDIR);
  pfo.set(fmh.get_sub(obdry_u),recv_o);
  
  // bulk
  for(int x=0;x<Nx-1;++x){
    for(int n=0;n<shiftSite_eo.slice_size(x,XDIR);++n)
      pfe.set(fmh.islice(shiftSite_eo.xslice(x,n,XDIR)),
	      fo[fmh.islice(shiftSite_oe.xslice(x+1,n,XDIR))]);
    for(int n=0;n<shiftSite_oe.slice_size(x,XDIR);++n)
      pfo.set(fmh.islice(shiftSite_oe.xslice(x,n,XDIR)),
	      fe[fmh.islice(shiftSite_eo.xslice(x+1,n,XDIR))]);
  }

  CCIO::cout<<"forward shifted fields:"<<std::endl;
  for(int i=0;i<ie.size();++i)
    CCIO::cout<<" (pfe["<<i<<"],pff["<<ie[i]<<"])=("<<pfe[i]<<","<<pff[ie[i]]<<")"
	      <<" (pfo["<<i<<"],pff["<<io[i]<<"])=("<<pfo[i]<<","<<pff[io[i]]<<")"
	      <<std::endl;  

  //////////////////////////////////////////////////////////////////
  // test of Dirac mult
  Field wf = D.mult(ff);
  Field we = Deo.mult_ee(fe); we+= Deo.mult_eo(fo);  
  Field wo = Deo.mult_oe(fe); wo+= Deo.mult_oo(fo);  

  for(int i=0;i<ie.size();++i)
    CCIO::cout<<" (we["<<i<<"],wf["<<ie[i]<<"])=("<<we[i]<<","<<wf[ie[i]]<<")"
	      <<" (wo["<<i<<"],wf["<<io[i]<<"])=("<<wo[i]<<","<<wf[io[i]]<<")"
	      <<std::endl;  
*/
  //------------ e/o preconditioned solver ------------------
  int    Niter= 2000;
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
  

  /*
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
  */


  return 0;
}

