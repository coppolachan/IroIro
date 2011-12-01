/*!
 * @file test_wilson_EvenOdd.cpp
 *
 * @brief Tests for the propagators 
 *
 */
#include "test_wilson_EvenOdd.hpp"
#include "Measurements/FermionicM/fermion_meas_factory.hpp"


#include "Dirac_ops/dirac_wilson_EvenOdd.h"
#include "Solver/solver_CG.h"
#include "Solver/solver_BiCGStab.h"
#include "Measurements/FermionicM/qprop_EvenOdd.h"
#include "Measurements/FermionicM/mesonCorrel.h"
#include "Tools/randNum_MT19937.h"
#include <stdio.h>

using namespace std;
using namespace Format;

int Test_Wilson_EvenOdd::run(){

  // Generating source vector 

  // local source
  vector<int> spos(4,0); 
  //vector<int> spos(4,1); 
  Source_local<Format_F> src(spos,CommonPrms::instance()->Nvol());

  /*
  // noise source
  unsigned long init[4]={0x123, 0x234, 0x345, 0x456};
  int length=4;
  RandNum_MT19937 rand(init, length);

  Source_wnoise<Format_F> src(rand,CommonPrms::instance()->Nvol());
  */

  /*
  //wall source not working
  Source_wall<Format_F> src(0,CommonPrms::instance()->Nvol());
  */

  prop_t sq; // propagator

  // Without factories -----------------------------------------------------
  //Dirac Kernel definition
  Dirac_Wilson         D(  1.0/6.0, &(conf_.U));
  Dirac_Wilson_EvenOdd Deo(1.0/6.0, &(conf_.U));

  SiteIndex_eo* ieo = SiteIndex_eo::instance();
  vector<int> esec= ieo->esec();
  vector<int> osec= ieo->osec();

  Format_F fmt(CommonPrms::instance()->Nvol());
  valarray<size_t> ie = fmt.get_sub(esec);
  valarray<size_t> io = fmt.get_sub(osec);

  // source generation
  Field ff= src.mksrc(1,1);
  Field fe= src.mksrc(esec,1,1);
  Field fo= src.mksrc(osec,1,1);
  /*
  // test of spritting into even/odd
  for(int i=0;i<ie.size();++i)
    CCIO::cout<<" (fe["<<i<<"],ff["<<ie[i]<<"])=("<<fe[i]<<","<<ff[ie[i]]<<")"
	      <<" (fo["<<i<<"],ff["<<io[i]<<"])=("<<fo[i]<<","<<ff[io[i]]<<")"
	      <<std::endl;  
  */

  /*
  // test of shift_field (up)
  Format_F fmh(CommonPrms::instance()->Nvol()/2);

  ShiftField_up<Format_F>      ffp(ff,&fmt,0);  valarray<double> ffpv(ffp.getva());
  ShiftField_even_up<Format_F> fep(fe,&fmh,0);  valarray<double> fepv(fep.getva());
  ShiftField_odd_up<Format_F>  fop(fo,&fmh,0);  valarray<double> fopv(fop.getva());

  for(int i=0;i<ie.size();++i)
    CCIO::cout<<" (fep["<<i<<"],ffp["<<io[i]<<"])=("<<fepv[i]<<","<<ffpv[io[i]]<<")"
	      <<" (fop["<<i<<"],ffp["<<ie[i]<<"])=("<<fopv[i]<<","<<ffpv[ie[i]]<<")"
	      <<std::endl;  

  // test of shift_field (dn)
  ShiftField_dn<Format_F>      ffm(ff,&fmt,0);  valarray<double> ffmv(ffm.getva());
  ShiftField_even_dn<Format_F> fem(fe,&fmh,0);  valarray<double> femv(fem.getva());
  ShiftField_odd_dn<Format_F>  fom(fo,&fmh,0);  valarray<double> fomv(fom.getva());

  for(int i=0;i<ie.size();++i)
    CCIO::cout<<" (fem["<<i<<"],ffm["<<io[i]<<"])=("<<femv[i]<<","<<ffmv[io[i]]<<")"
	      <<" (fom["<<i<<"],ffm["<<ie[i]<<"])=("<<fomv[i]<<","<<ffmv[ie[i]]<<")"
	      <<std::endl;  
  */
  /*
  // test of Dirac mult
  Field wf = D.mult(ff);
  Field we = Deo.mult_ee(fe); we+= Deo.mult_eo(fo);  
  Field wo = Deo.mult_oe(fe); wo+= Deo.mult_oo(fo);  

  for(int i=0;i<ie.size();++i)
    CCIO::cout<<" (we["<<i<<"],wf["<<ie[i]<<"])=("<<we[i]<<","<<wf[ie[i]]<<")"
	      <<" (wo["<<i<<"],wf["<<io[i]<<"])=("<<wo[i]<<","<<wf[io[i]]<<")"
	      <<std::endl;  
  */
  
  // test of Solver & Qprop
  int    Niter= 1000;
  double stop_cond = 1.0e-24;
  Fopr_DdagD DdagD(&Deo);
  //Solver_BiCGStab solver(stop_cond,Niter,&DdagD);
  Solver_CG solver(stop_cond,Niter,&DdagD);  

  Qprop_EvenOdd qprop(&Deo,&solver);
  qprop.calc(sq,src);


  //  Using factories ------------------------------------------------------
  /*
  QuarkPropagator* QP;
  XML::descend(Wilson_EO_node_, "QuarkPropagator");

  QuarkPropagatorFactory* QP_Factory = 
    QuarkPropagators::createQuarkPropagatorFactory(Wilson_EO_node_);

  QP = QP_Factory->getQuarkProp(conf_);
  QP->calc(sq,src);
  */
  //------------------------------------------------------------------------

  CCIO::cout<<"quark propagator obtained"<<std::endl;
  
  // meson correlators
  MesonCorrel<Format_F> meson;
  vector<double> mcorr = meson.pp(sq,sq);  
  vector<double>::const_iterator it=mcorr.begin();
  int t=0;
  while(it!=mcorr.end()) pprintf ("%d %.8e\n",t++, *it++);
  return 0;
}

