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
#include "Solver/solver_CG.h"
//#include "Solver/solver_BiCGStab.h"
#include "Measurements/FermionicM/qprop_optimalDomainWall.hpp"
#include "Measurements/FermionicM/mesonCorrel.h"
#include "Measurements/FermionicM/source.h"

#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <ctime>

using namespace std;
using namespace Format;

const Field Test_ResMass::delta(const Dirac_optimalDomainWall_4D& DWF, Field& phi){
  //Delta function = 1/4 * (1- sign^2(Hw))
  Field sign = DWF.signKernel(phi);
  Field delta = DWF.signKernel(sign); //sign^2(Hw)
  delta -= phi;  //sign^2(Hw) -1
  
  delta *= -0.25; // 1/4*(1-sign^2(Hw))
  
  return delta;
}

int Test_ResMass::run(XML::node node) {
  std::ifstream input_file;
  input_file.open("param.dat");

  XML::descend(node, "DiracOperator");
  // operator
  double mzero;
  int N5d;    // the length in the 5th direction (must be even)
  double b;   // scale factor (b!=1 for scaled Shamir H_T)
  double c;   // the kernel (H_W (c=0) or H_T (c=1))
  /* Kernel is given by H=gamma_5 b D_W /( 2 + c D_W ) */
  int approx; // the approx (tanh (approx=0) or Zolo (approx=1))
  double mq;
  input_file >> mzero;
  input_file >> N5d;
  input_file >> b;
  input_file >> c;
  input_file >> approx;
  input_file >> mq;

  CCIO::cout << "Residual mass testing with" << endl;
  CCIO::cout << " mzero  = " << mzero << endl;
  CCIO::cout << " N5d    = " << N5d << endl;
  CCIO::cout << " b      = " << b << " (scale factor)" << endl;
  CCIO::cout << " c      = " << c << " (H_W (c=0) or H_T (c=1))" << endl;
  CCIO::cout << " approx = " << approx << " (tanh (approx=0) or Zolo (approx=1))" << endl;
  CCIO::cout << " mq     = " << mq << endl;

  Dirac_Wilson* Kernel = new Dirac_Wilson(mzero, &(conf_.U));

  vector<double> omega(N5d);
  if (approx == 0) {
    for (int s = 0; s < N5d; ++s) omega[s] = 1.0;
  } else {
    double lambda_min, lambda_max;
    input_file >> lambda_min;
    input_file >> lambda_max;
    CCIO::cout << " lambda_min = " << lambda_min << endl;
    CCIO::cout << " lambda_max = " << lambda_max << endl;
    omega = DomainWallFermions::getOmega(N5d,lambda_min,lambda_max);
  }

  CCIO::cout << "Creating the DWF operator" << endl;
  Dirac_optimalDomainWall Ddwf_5d(b,c,mq,omega,Kernel);
  
  // quark propagator
  double stop_cond = 1.0e-10;
  int Niter= 20000;
  CCIO::cout << " stop_cond = " << stop_cond << endl;
  CCIO::cout << " Niter     = " << Niter << endl;

  // standard construction (factories will use a similar one)
  Dirac_optimalDomainWall Ddwf_PV(Ddwf_5d, PauliVillars);
  Solver* SolvDWF = new Solver_CG(stop_cond,Niter,new Fopr_DdagD(&Ddwf_5d));
  Solver* SolvPV  = new Solver_CG(stop_cond,Niter,new Fopr_DdagD(&Ddwf_PV));
  //  Solver* SolvDWF 
  //   = new Solver_BiCGStab(stop_cond,Niter,new Fopr_DdagD(&Ddwf_5d));
  //  Solver* SolvPV  
  //   = new Solver_BiCGStab(stop_cond,Niter,new Fopr_DdagD(&Ddwf_PV));
  Dirac_optimalDomainWall_4D DiracDWF_4d(Ddwf_5d,SolvDWF,SolvPV);
  QpropDWF QuarkPropagator(DiracDWF_4d);

  vector<int> spos(4,0); 
  //Source generator
  Source_local<Format_F> src(spos,CommonPrms::instance()->Nvol());

  prop_t sq;  //Defines a vector of fields
  CCIO::cout << "Calculating propagator\n";
  QuarkPropagator.calc(sq,src);
 
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
      CCIO::cout << "Numerator = ("<<mres_numerator<<","<<im_check<<")\n";
      
      //Denominator
      Denom = sq[c+3*s];
      Denom -= src.mksrc(s,c); // (D^-1 - 1)*src
      Denom /= (1.0 - mq);
      mres_denominator += Denom*Denom;
      CCIO::cout << "Denominator = " << mres_denominator << endl;
      CCIO::cout << "Residual mass = " << mres_numerator/mres_denominator << endl;
    }
  }

  return 0;
}
