/*!
 * @file test_MultiShiftSolver.cpp
 *
 * @brief Declaration of functions to test MultiShiftSolver_CG 
 */
#include "test_MultiShiftSolver.hpp"
#include "Solver/multiShiftSolver_CG.hpp"
#include "include/format_F.h"
#include "include/fopr.h"
#include "Measurements/FermionicM/source_types.hpp"
#include "Measurements/FermionicM/qprop_MultiShift.hpp"
#include "Dirac_ops/dirac_wilson.hpp"
#include "Dirac_ops/dirac_clover.hpp"
#include "Dirac_ops/dirac_DomainWall.hpp"


#include <vector>
#include <stdio.h>
#include <iostream>

using namespace std;

int Test_MultiShiftSolver::run(){
  test1();
  return 0;
}


int Test_MultiShiftSolver::test1(){

  // Definition of shifts
  prop_t  xqs;
  vector<double> mass_shifts;
  mass_shifts.push_back(0.05);
  mass_shifts.push_back(0.10);
  mass_shifts.push_back(0.20);
  mass_shifts.push_back(0.30);

  int N5d   = 6;
  double M0 = -1.6;
  double c  = 1.0;
  double b  = 1.0;
  double mq = 0.02;
  vector<double> omega(N5d,1.0);

  Dirac_optimalDomainWall* Kernel = new Dirac_optimalDomainWall(b,c,M0,mq,omega,&(Gauge.data));


  // Definition of source 
  vector<int> spos(4,0);
  Source_local<Format::Format_F> Source(spos,
					CommonPrms::instance()->Nvol()*N5d);

  // Definition of Dirac Kernel
  //Dirac* Kernel = new Dirac_Wilson(M0, &(Gauge.data));
  //Dirac* Kernel = new Dirac_Clover(0.01, 1.0, &(Gauge.data));


  // Definition of the Solver
  int    Niter= 1000;
  double stop_cond = 1.0e-24;
  MultiShiftSolver* Solver = 
    new MultiShiftSolver_CG(new Fopr_DDdag(Kernel),
			    stop_cond,
			    Niter);
  

  // Solver test
  xqs.resize(mass_shifts.size());
  for (int i = 0; i < mass_shifts.size(); ++i)
    xqs.push_back(Field(CommonPrms::instance()->Nvol()*N5d));
  double residual;
  int Nconv;
  Solver->solve(xqs, Source.mksrc(0,0), mass_shifts, residual, Nconv);



  
  // Quark Propagator test
  /*
    Qprop_MultiShift QuarkPropagator(Kernel, Solver);
    QuarkPropagator.calc(xqs,
    Source,
    mass_shifts);
  */
  return 0;
}
