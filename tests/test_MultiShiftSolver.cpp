/*!
 * @file test_MultiShiftSolver.cpp
 *
 * @brief Declaration of functions to test MultiShiftSolver_CG 
 */
#include "test_MultiShiftSolver.h"
#include "Solver/multiShiftSolver_CG.hpp"
#include "include/format_F.h"
#include "include/fopr.h"
#include "Measurements/FermionicM/source_types.hpp"
#include "Measurements/FermionicM/qprop_MultiShift.hpp"
#include "Dirac_ops/dirac_wilson.h"

#include <vector>
#include <stdio.h>
#include <iostream>

using namespace std;

int Test_MultiShiftSolver::run(){
  test1();
  return 0;
}


int Test_MultiShiftSolver::test1(){

  Format::Format_F ff(CommonPrms::instance()->Nvol());

  // Definition of shifts
  prop_t  xqs;
  vector<double> mass_shifts;
  mass_shifts.push_back(0.10);
  mass_shifts.push_back(0.20);
  mass_shifts.push_back(0.30);
  mass_shifts.push_back(0.40);

  // Definition of source 
  vector<int> spos(4,0);
  Source_local<Format::Format_F> Source(spos,
					CommonPrms::instance()->Nvol());

  // Definition of Dirac Kernel
  double M0=1.6;
  Dirac* Kernel = new Dirac_Wilson(-M0, &(Gauge.U));

  // Definition of the Solver
  int    Niter= 1000;
  double stop_cond = 1.0e-24;
  MultiShiftSolver* Solver = 
    new MultiShiftSolver_CG(new Fopr_DDdag(Kernel),
			    stop_cond,
			    Niter);
  
  
  Qprop_MultiShift QuarkPropagator(Kernel, Solver);
  QuarkPropagator.calc(xqs,
		       Source,
		       mass_shifts);

  return 0;
}
