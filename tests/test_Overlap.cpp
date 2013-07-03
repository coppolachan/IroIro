/*!
 * @file test_Overlap.hpp
 *
 * @brief Definition of functions to test the overlap functions
 */
#include "test_Overlap.hpp"
#include "Dirac_ops/dirac_wilson.hpp"
#include "Dirac_ops/dirac_overlap_Zolotarev.hpp"
#include "EigenModes/eigenProc_Zolotarev.h"
#include "Solver/solver_CG.hpp"
#include "Measurements/FermionicM/qprop.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include "Measurements/FermionicM/source_types.hpp"
#include "Fields/field_expressions.hpp"

#include <stdio.h>


int Test_Overlap::run() {
  sign_subt();
  ovsolver_subt();
}

int Test_Overlap::sign_subt(){
  using namespace FieldExpression;
  std::cout << "Test_Overlap::sign_subt" << std::endl;

  // test of sign function
  std::vector<int> spos(4,0); 
  Source_local<Format::Format_F> source(spos,CommonPrms::instance()->Nvol());

  // Creation of Kernel Dw
  double M0 = 1.6;
  DiracWilsonLike* Kernel = new Dirac_Wilson(-M0,&(Gauge.data));

  EigenPrms Eprms(//                     (calc of eigenmodes)
		  100,                   // Nmm
		  20,20,                 // Nk_low, Nk_high 
		  50,50,                 // Np_low, Np_high 
		  500,500,               // Niter_low, Niter_high
		  1.0e-22,1.0e-20,       // enorm_low, enorm_high
		  0.15,10.0,             // vthrs_low, vthrs_high

		  //                     (Zolotarev approx.)    
		  1.0e-20,               // stp_cnd
		  1000,                  // Niter
		  16);                   // Npoly  
  
  EigenProc_Zolotarev eproc(Kernel,Eprms);

  EigenData ed;
  eproc.calc(ed);
  
  Fopr_signH_Zolotarev signHw(Kernel,Eprms,&ed);
  
  Field xq(signHw.fsize());
    for(int s = 0; s < CommonPrms::instance()->Nd(); ++s){
      for(int c = 0; c < CommonPrms::instance()->Nc(); ++c){
      Field b = source.mksrc(s,c);
      xq= signHw.mult_dag(signHw.mult(b)) -b;
      printf(" |sign^2(b) - b| = %16.8e\n",xq.norm());
    }
  }
  return 0;
}

int Test_Overlap::ovsolver_subt(){
  using namespace FieldExpression;

  // Definition of Kernel
  double M0 = 1.6;
  double mq = 0.15;
  DiracWilsonLike* Kernel = new Dirac_Wilson(-M0,&(Gauge.data));

  EigenPrms  Eprms(//                     (calc of eigenmodes)
		   100,                   // Nmm
		   20,20,                 // Nk_low, Nk_high 
		   50,50,                 // Np_low, Np_high 
		   500,500,               // Niter_low, Niter_high
		   1.0e-22,1.0e-20,       // enorm_low, enorm_high
		   0.15,10.0,             // vthrs_low, vthrs_high
		   //                     (Zolotarev approx.)    
		   1.0e-20,               // stp_cnd
		   1000,                  // Niter
		   16);                   // Npoly  
  
  // construction of Overlap operator
  EigenProc_Zolotarev eproc(Kernel,Eprms);

  EigenData ed;
  eproc.calc(ed);

  Fopr_signH_Zolotarev signHw(Kernel,Eprms,&ed);  
  Dirac* Overlap = new Dirac_overlap_Zolotarev(M0,mq,&signHw);
  
  std::cout << "Dirac_overlap_Zolotarev  was constructed\n";

  // construction of Qprop
  Solver* Solv = new Solver_CG(Eprms.stp_cnd, 
			       Eprms.Niter, 
			       new Fopr_DdagD(Overlap));
  Qprop qprop(Overlap, Solv);
  std::cout << "Qprop was constructed\n";

  // source
  std::vector<int> spos(4,0); 
  Source_local<Format::Format_F> src(spos,CommonPrms::instance()->Nvol());

  // getting quark propagator
  std::vector<Field> sq;
  qprop.calc(sq,src);
  
  GammaMatrices::Unit Gamma;//pion
  MesonCorrelator meson(Pion);
  std::vector<double> mcorr = meson.calculate< Format::Format_F >(sq,sq);
  for(int t = 0; t < mcorr.size(); ++t)
    std::cout << t << "  "<< mcorr[t] << std::endl;
  
  return 0;
}

