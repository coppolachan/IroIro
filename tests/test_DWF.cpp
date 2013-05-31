/*!
 * @file test_DWF.cpp
 * @brief Tests for the propagators 
 */
#include <stdio.h>

#include "test_DWF.hpp"
#include "Dirac_ops/dirac_wilson_Brillouin.hpp"
#include "Dirac_ops/dirac_DomainWall.hpp"
#include "Dirac_ops/dirac_DomainWall_EvenOdd.hpp"
#include "Dirac_ops/dirac_DomainWall_4D_eoSolv.hpp"
#include "Dirac_ops/dirac_DomainWall_4D_fullSolv.hpp"
#include "Solver/solver_CG.hpp"
#include "Solver/solver_BiCGStab.hpp"
#include "Measurements/FermionicM/qprop_DomainWall.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include "Tools/RandomNumGen/randNum_MT19937.h"
#include "Measurements/FermionicM/source_types.hpp"
#include "Measurements/GaugeM/staples.hpp"

using namespace std;
using namespace Format;
using namespace DomainWallFermions;
using namespace EvenOddUtils;

int Test_DWF::run(){
  CCIO::cout<<"Test_DWF::run() called\n";

  Staples stpl;
  double plq = stpl.plaquette(Gfield_);
  CCIO::cout<<" plaq="<<plq<<std::endl;

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
  // source generation
  Field ff= src.mksrc(1,1);

  /************************************************************************************/
  Field* u = &(Gfield_.data);
  double M0 = -1.6;

  // creation of Dirac_Wilson operators 
  Dirac_Wilson Dw(   M0,u);
  Dirac_Wilson Dw_eo(M0,u,Dop::EOtag());
  Dirac_Wilson Dw_oe(M0,u,Dop::OEtag());
  /*
  Dirac_Wilson_Brillouin Db(M0,u);
  CCIO::cout<<"Dbrillouin\n";
  Field tmp = Db.mult(ff);
  */
  int N5=4;
  double b=2.0;
  double c=0.0;
  double mq =0.05;
  ffmt_t fmt5d(CommonPrms::instance()->Nvol(),N5);

  std::vector<double> omega(N5,1.0);

  double prec = 1.0e-16;
  int max_iter = 600;

  ///////
  Dirac_optimalDomainWall Ddwf(b,c,M0, mq,omega,&Dw,u);
  Dirac_optimalDomainWall Ddpv(b,c,M0,1.0,omega,&Dw,u);

  Fopr_DdagD DdagDdwf(&Ddwf);
  Fopr_DdagD DdagDdpv(&Ddpv);

  Solver_CG slv_dwf(prec,max_iter,&DdagDdwf);
  Solver_CG slv_dpv(prec,max_iter,&DdagDdpv);

  Dirac_optimalDomainWall_4D_fullSolv D4f( &Ddwf,&Ddpv,&slv_dwf,&slv_dpv);
  //Dirac_optimalDomainWall_4D_fullSolv D4f(&Ddwf,prec,prec,max_iter);

  Dirac_optimalDomainWall_EvenOdd Ddwf_eo(b,c,M0, mq,omega,&Dw_eo,&Dw_oe,u);
  Dirac_optimalDomainWall_EvenOdd Ddpv_eo(b,c,M0,1.0,omega,&Dw_eo,&Dw_oe,u);

  Fopr_DdagD DdagDdwf_eo(&Ddwf_eo);
  Fopr_DdagD DdagDdpv_eo(&Ddpv_eo);

  Solver_CG slv_dwf_eo(prec,max_iter,&DdagDdwf_eo);
  Solver_CG slv_dpv_eo(prec,max_iter,&DdagDdpv_eo);

  Inverter_WilsonLike invDdwf(&Ddwf_eo,&slv_dwf_eo);
  Inverter_WilsonLike invDdpv(&Ddpv_eo,&slv_dpv_eo);

  Dirac_optimalDomainWall_4D_eoSolv D4eo(N5,mq,&invDdwf,&invDdpv);

  /*--------mult test---------*/

  
 // plain Wilson fermion
  Field Sdw(Dw.fsize()); Sdw.set(0,1.0);
  Field Wwf =  Dw.mult(Sdw);
  double Nwf = Wwf.norm();
  CCIO::cout<<"Nwf="<<Nwf<<"\n";

  // e/o Wilson fermion
  Field Sweo(Dw_eo.fsize()); Sweo.set(0,1.0);
  Field Wweo =  Dw_eo.mult(Sweo);
  double Nweo = Wweo.norm();
  CCIO::cout<<"Nweo="<<Nweo<<"\n";

 // plain 5D DWF
  Field Sdwf(Ddwf.fsize()); 
  Sdwf.set(fmt5d.index( 0,0,0),0.5);
  Sdwf.set(fmt5d.index( 0,0,3),0.5);
  Sdwf.set(fmt5d.index(12,0,3),0.5);

  Field Wdwf = Ddwf.mult(Sdwf);
  double Ndwf =Wdwf.norm();
  Wdwf = Ddwf.mult_dag(Sdwf);
  double Ndwfd =Wdwf.norm();
  CCIO::cout<<"Ndwf="<<Ndwf<<" Ndwfd="<<Ndwfd<<"\n";

 // e/o 5D DWF
  Field Sdweo(Ddwf_eo.fsize()); Sdweo.set(0,1.0);
  Field Wdweo = Ddwf_eo.mult(Sdweo);
  double Ndweo =Wdweo.norm();
  CCIO::cout<<"Ndweo="<<Ndweo<<"\n";

  double nwfl=0.0; 
  double nweo=0.0; 

  CCIO::cout<<"test of 4D full solver \n";
  Field S4(D4f.fsize()); S4.set(0,1.0);
  //Field wfl = D4f.mult(ff);
  Field wfl = D4f.mult(S4);
  nwfl = wfl.norm();

  CCIO::cout<<"test of 4D e/o solver \n";
  Field S4eo(D4eo.fsize()); S4eo.set(0,1.0);  
  //  Field weo = D4eo.mult(ff); 
  Field weo = D4eo.mult(S4eo); 
  nweo = weo.norm();

  CCIO::cout<<"nwfl="<<nwfl<<" nweo="<<nweo<<"\n";

  //////////////////////////////////////////////////////

  SiteIndex_EvenOdd* ieo = SiteIndex_EvenOdd::instance();
  vector<int> esec= ieo->esec();
  vector<int> osec= ieo->osec();

  Format::Format_F fmt(CommonPrms::instance()->Nvol());
  valarray<size_t> ie = fmt.get_sub(esec);
  valarray<size_t> io = fmt.get_sub(osec);

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

