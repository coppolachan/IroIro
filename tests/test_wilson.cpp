/*!
 * @file test_wilson.cpp
 * @brief Tests for the propagators 
 */
#include "test_wilson.hpp"
#include "Dirac_ops/dirac_wilson.hpp"
#include "Solver/solver_CG.hpp"
#include "Solver/solver_BiCGStab.hpp"
#include "Measurements/FermionicM/quark_prop_meas_factory.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include "Measurements/FermionicM/qprop.hpp"
#include "Measurements/FermionicM/source_types.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include "Tools/RandomNumGen/randNum_MT19937.h"
#include <stdio.h>
#include <time.h>
using namespace std;
using namespace Format;

timespec diff(timespec start, timespec end)
{
	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}

int Test_Wilson::run(){
  prop_t sq;// propagator
  // Use factories to construct the propagator
  QuarkPropagator* QP;
  XML::descend(wilson_node_, "QuarkPropagator");

  Staples Staple;
  CCIO::cout << "Plaquette : " << Staple.plaquette(conf_) << std::endl;

   QuarkPropagatorFactory* QP_Factory = 
    QuarkPropagators::createQuarkPropagatorFactory(wilson_node_);
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
  Dirac_Wilson* Kernel = new Dirac_Wilson(-1.15, &(conf_.data));
  //Dirac* Kernel = new Dirac_Clover(1.0/6.0, 1.0, &(conf_.U));



  // Test for mult and mult_new

  int calls = 100;
  timespec ts1, ts2;

  Field const_src(Kernel->fsize(),1.0);
  const_src/= const_src.norm();

  Field v2(Kernel->fsize());
  clock_gettime(CLOCK_REALTIME, &ts1);
  for (int i = 0; i < calls; ++i) {
    //    v2 = Kernel->mult(src.mksrc(0,0));
    v2 = Kernel->mult(const_src);
    double v2nrm = v2.norm();
    CCIO::cout<<"v2norm="<<v2nrm<<"\n";
  }
  clock_gettime(CLOCK_REALTIME, &ts2);
  cout<<diff(ts1,ts2).tv_sec<<"."<<diff(ts1,ts2).tv_nsec<<endl;
  /*
  clock_gettime(CLOCK_REALTIME, &ts1);
  for (int i = 0; i < calls; ++i) {
    Kernel->mult_new(v2,src.mksrc(0,0));
  }
  clock_gettime(CLOCK_REALTIME, &ts2);
  cout<<diff(ts1,ts2).tv_sec<<"."<<diff(ts1,ts2).tv_nsec<<endl;
  */



  /*  
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
  CCIO::cout << "sq size " << sq.size()<< "\n";
  for (int i = 0; i < sq.size(); i++){
    Field temp = sq[i];
    for (int j = 0; j < temp.size(); j++){
      CCIO::cout << "["<<j<<"]["<<i<<"] "<<temp[j] <<  "\n";
    }  
  }

  //  QP = QP_Factory->getQuarkProp(conf_);
  //  QP->calc(sq,src);
  
  CCIO::cout<<"quark propagator obtained"<<std::endl;
  */
  return 0;
}

