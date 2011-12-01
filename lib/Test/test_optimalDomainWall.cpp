/*!
 * @file test_optimalDomainWall.cpp
 *
 * @brief Definition of classes for testing the Dirac_optimalDomainWall classes
 *
 */
#include "test_optimalDomainWall.hpp"
#include "Communicator/comm_io.hpp"
#include "Measurements/FermionicM/fermion_meas_factory.hpp"
#include "include/fopr.h"
#include "Tools/randNum_MT19937.h"
//#include "Solver/solver_CG.h"
#include "Solver/solver_BiCGStab.h"
#include "Measurements/FermionicM/qprop_optimalDomainWall.hpp"
#include "Measurements/FermionicM/mesonCorrel.h"

#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <ctime>

using namespace std;
using namespace Format;


int Test_optimalDomainWall::mult5d_test(Dirac_optimalDomainWall& DWF5d,
					Field& InputField,
					int iterations){
  Field dphi;
  for (int i = 0; i < iterations; i++)
    dphi = DWF5d.mult(InputField);
  double dphi_norm = dphi.norm();    // dphi= Ddwf*phi
  CCIO::cout << "dphi.norm: " << dphi_norm << endl;   
}

int Test_optimalDomainWall::mult5d_dag_test(Dirac_optimalDomainWall& DWF5d,
					    Field& InputField,
					    int iterations){
  Field dphi;
  for (int i = 0; i < iterations; i++)
    dphi = DWF5d.mult_dag(InputField);
  double dphi_norm = dphi.norm();    // dphi= Ddwf*phi
  CCIO::cout << "dphi.norm: " << dphi_norm << endl;   
}

int Test_optimalDomainWall::mult5d_gamma5_test(Dirac_optimalDomainWall& DWF5d,
					       Field& InputField,
					       int iterations){
  Field dphi;
  Field g5psi,dg5psi,g5dg5psi;
  for (int i = 0; i < iterations; i++){
    g5psi = DWF5d.R5g5(InputField);         
    dg5psi = DWF5d.mult(g5psi);
    g5dg5psi = DWF5d.R5g5(dg5psi);
  }

  double g5dg5psi_norm = g5dg5psi.norm();

  Field vdiff = g5dg5psi;
  Field ddagphi = DWF5d.mult_dag(InputField);
  vdiff -= ddagphi;
  double vdiff_norm = vdiff.norm();

  CCIO::cout << "(g5 d g5 v) norm: " << g5dg5psi_norm << endl;
  CCIO::cout << "dhc v  norm: " << ddagphi.norm() << endl;
  CCIO::cout << "Difference norm: " << vdiff_norm << endl;
}

int Test_optimalDomainWall::run(XML::node node){
  Dirac_optimalDomainWall* DiracODWF;

  // operator without factories
  int N5d = 6;
  double mzero = -1.8;
  double c = 0.0;
  double b = 2.0;
  double mq = 0.0;
  vector<double> omega(N5d,1.0);
  Dirac_Wilson* Kernel = new Dirac_Wilson(mzero, &(conf_.U));
  Dirac_optimalDomainWall Ddwf_5d(b, c, mq, omega, Kernel);
  /////////////////////////////

  // operator using factories
  XML::descend(node, "DomainWall");
  DiracDWF5dFactory* DWF_Factory = new DiracDWF5dFactory(node);
  DiracODWF = DWF_Factory->getDiracOperator(&(conf_.U));
  ///////////////////////////////////////

  // prepare pf field
  unsigned long init[4]={0x123, 0x234, 0x345, 0x456};
  int length=4;
  RandNum_MT19937 rand(init,length);

  valarray<double> vphi(Ddwf_5d.fsize());
  rand.get(vphi);
  Field phi(vphi);    // phi: generated from random numbers
  double phi_norm = phi.norm();
  CCIO::cout << "phi,size,norm: " << phi.size() << " " << phi_norm << endl;
  
  valarray<double> vpsi(Ddwf_5d.fsize());
  rand.get(vpsi);
  Field psi(vpsi);
  double psi_norm = psi.norm(); // psi: generated from randum numbers
  CCIO::cout << "psi,size,norm: " << psi.size() << " " << psi_norm << endl;

  // mult test
  clock_t start;
  CCIO::cout << ".::: Test Dirac_optimalDomainWall.mult(f)\n";
  start = clock();
  mult5d_test(*DiracODWF,phi,100);
  CCIO::cout<< "Time for 100 calls : " 
      << ( ( clock() - start ) / (double)CLOCKS_PER_SEC ) <<'\n';

  CCIO::cout << ".::: Test Dirac_optimalDomainWall.mult_dag(f)\n";
  start = clock();
  mult5d_dag_test(*DiracODWF,phi,100);
  CCIO::cout<< "Time for 100 calls : " 
      << ( ( clock() - start ) / (double)CLOCKS_PER_SEC ) <<'\n';

  // Gamma 5 Hermiticity
  CCIO::cout << ".::: Test Dirac_optimalDomainWall gamma5 hermiticity " 
       << " vdiff = ( G5 D G5 - D.hc ) psi \n";
  start = clock();
  mult5d_gamma5_test(*DiracODWF,psi,100);
  CCIO::cout<< "Time for 100 calls : " 
      << ( ( clock() - start ) / (double)CLOCKS_PER_SEC ) <<'\n';
  
  // quark propagator
  double stop_cond = 1.0e-24;
  int    Niter= 1000;

  /* It follows a standard construction (factories will use a similar one)
  Dirac_optimalDomainWall Ddwf_PV(Ddwf_5d, PauliVillars);
  Solver* SolvDWF = new Solver_BiCGStab(stop_cond,Niter,new Fopr_DdagD(&Ddwf_5d));
  Solver* SolvPV = new Solver_BiCGStab(stop_cond,Niter,
				       new Fopr_DdagD(&Ddwf_PV));
  Dirac_optimalDomainWall_4D DiracDWF_4d(Ddwf_5d,*SolvDWF, *SolvPV);
  QpropDWF QuarkPropagator(DiracDWF_4d);
  //////////////////////////////////// */

  // Here uses the default constructor with default solver
  QpropDWF QuarkPropagator(Ddwf_5d,
   			   stop_cond,
   			   Niter);
  
  vector<int> spos(4,0); 
  Source_local<Format_F> src(spos,CommonPrms::instance()->Nvol());

  prop_t sq;
  QuarkPropagator.calc(sq,src);

  // conf_.U.set(0,1.0);//change configuration - test
  // QuarkPropagator.calc(sq,src,0,0);

  return 0;
}
