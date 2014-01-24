/*!
 * @file test_ResidualMass.cpp
 *
 * @brief Definition of classes for calculating the Residual Mass
 *
 */
#include "test_ResidualMass.hpp"
#include "Communicator/comm_io.hpp"
#include "Communicator/fields_io.hpp"

#include "Tools/randNum_MT19937.h"
//#include "Solver/solver_CG.h"
#include "Solver/solver_BiCGStab.h"
#include "Measurements/FermionicM/qprop_DomainWall.hpp"
#include "Measurements/FermionicM/mesonCorrel.h"
#include "Measurements/FermionicM/source.h"

#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <ctime>

using namespace std;
using namespace Format;

const Field Test_ResMass::delta(const Dirac_DomainWall_4D& DWF, Field& phi){
  //Delta function = 1/4 * (1- sign^2(Hw))
  Field sign = DWF.signKernel(phi);
  Field delta = DWF.signKernel(sign); //sign^2(Hw)
  delta -= phi;  //sign^2(Hw) -1
  
  delta *= -0.25; // 1/4*(1-sign^2(Hw))
  
  return delta;
}

int Test_ResMass::run(XML::node node) {
  // operator
  int N5d = 6;
  double mzero = -1.8;
  double c = 0.0;
  double mq = 0.0;
  vector<double> omega(N5d,1.0);

  Dirac_Wilson* Kernel = new Dirac_Wilson(mzero, &(conf_.U));
 
  Dirac_DomainWall Ddwf_5d(c, mq, omega, Kernel);
  
  // quark propagator
  double stop_cond = 1.0e-24;
  int    Niter= 1000;

  //It follows a standard construction (factories will use a similar one)
  Dirac_DomainWall Ddwf_PV(Ddwf_5d, PauliVillars);
  Solver* SolvDWF = new Solver_BiCGStab(stop_cond,Niter,new Fopr_DdagD(&Ddwf_5d));
  Solver* SolvPV  = new Solver_BiCGStab(stop_cond,Niter,new Fopr_DdagD(&Ddwf_PV));
  Dirac_DomainWall_4D DiracDWF_4d(Ddwf_5d,SolvDWF,SolvPV);
  QpropDWF QuarkPropagator(DiracDWF_4d);
  //////////////////////////////////// 

  vector<int> spos(4,0); 
  //Source generator
  Source_local<Format_F> src(spos,CommonPrms::instance()->Nvol());

  prop_t sq;  //Defines a vector of fields
  CCIO::cout << "Calculating propagator\n";
  QuarkPropagator.calc(sq,src);
  //save the quark propagator somewhere
  CCIO::cout << "Saving propagator\n";
  CCIO::SaveOnDisk(sq, "./test_propagator.prop");

  // Cycle among Dirac and color indexes and contract
  // D^-1 * Delta * D^-1
  double mres_numerator = 0;
  double im_check = 0;
  double mres_denominator = 0;
  Field Delta, Denom;

  for (int s = 0; s < 4; s++) {
    for (int c = 0; c < 3; c++) {
      Delta = delta(DiracDWF_4d,sq[c+3*s]); // (Delta * D^-1)*source
      //Contracts
      mres_numerator += sq[c+3*s]*Delta;          // Re(sq[],Delta)    sq[]=D^-1*source
      im_check       += sq[c+3*s].im_prod(Delta); //should be always zero (just a check)
      CCIO::cout << "Contraction = ("<<mres_numerator<<","<<im_check<<")\n";
      
      //Denominator
      Denom = sq[c+3*s];
      Denom -= src.mksrc(s,c); // (D^-1 - 1)*src
      Denom /= (1.0 - mq);
      mres_denominator += Denom*Denom;
      
      CCIO::cout << "Residual mass = " << mres_numerator/mres_denominator << endl;
    }
  }


  return 0;
}
