/*!
 * @file test_optimalDomainWall.cpp
 *
 * @brief Definition of classes for testing the Dirac_optimalDomainWall classes
 *
 */
#include "test_DomainWall.hpp"
#include "Communicator/comm_io.hpp"
#include "Measurements/FermionicM/fermion_meas_factory.hpp"
#include "include/fopr.h"
#include "Tools/randNum_MT19937.h"
//#include "Solver/solver_CG.h"
#include "Solver/solver_BiCGStab.h"
#include "Measurements/FermionicM/qprop_DomainWall.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include "Measurements/FermionicM/source_types.hpp"

#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <ctime>

using namespace std;
using namespace Format;
#include <typeinfo> 

int Test_optimalDomainWall::
mult5d_test(const Dirac_optimalDomainWall& DWF5d,
	    const Field& phi,int times){

  Field Dphi(phi.size());
  for(int i=0; i<times; ++i) Dphi = DWF5d.mult(phi);

  double Dphi_norm = Dphi.norm();    // Dphi= Ddwf*phi
  CCIO::cout << "Dphi.norm: " << Dphi_norm << endl;   
}

int Test_optimalDomainWall::
mult5d_dag_test(const Dirac_optimalDomainWall& DWF5d,
		const Field& phi,int times){
  
  Field Dphi(phi.size());
  for(int i=0; i<times; ++i) Dphi = DWF5d.mult_dag(phi);

  double Dphi_norm = Dphi.norm();    // Dphi= Ddwf*phi
  CCIO::cout << "Dphi.norm: " << Dphi_norm << endl;   
}

int Test_optimalDomainWall::
mult5d_gamma5_test(const Dirac_optimalDomainWall& DWF5d,
		   const Field& phi,int times){

  Field g5psi(phi.size());
  Field Dg5psi(phi.size());
  Field g5Dg5psi(phi.size());

  for(int i=0; i<times; ++i){
    g5psi = DWF5d.R5g5(phi);         
    Dg5psi = DWF5d.mult(g5psi);
    g5Dg5psi = DWF5d.R5g5(Dg5psi);
  }

  double g5Dg5psi_norm = g5Dg5psi.norm();

  Field vdiff = g5Dg5psi;
  Field Ddagphi = DWF5d.mult_dag(phi);
  vdiff -= Ddagphi;
  double vdiff_norm = vdiff.norm();

  CCIO::cout << "||g5*D*g5*v||: " << g5Dg5psi_norm << std::endl;
  CCIO::cout << "||Ddag*v||: "    << Ddagphi.norm() << std::endl;
  CCIO::cout << "||Difference||: "<< vdiff_norm << std::endl;
}

int Test_optimalDomainWall::run(XML::node node){
  Dirac_optimalDomainWall* DiracODWF;

  // operator without factories
  int N5d = 8;
  double M0 = -1.6;
  double c = 1.0;
  double b = 1.0;
  double mq = 0.02;
  vector<double> omega(N5d,1.0);

  Dirac_optimalDomainWall Ddwf_5d(b,c,M0,mq,omega,&(conf_.U));
  /////////////////////////////
  XML::node QuarkProp_node = node;
  // operator using factories
  XML::descend(node, "DomainWall");
  DiracDWF5dFactory DWF_Factory(node);
  DiracODWF = DWF_Factory.getDiracOperator(&(conf_.U));

  
  XML::descend(QuarkProp_node, "QuarkPropagator");
  QPropDWFFactory  QP_DomainWallFact(QuarkProp_node);//uses specific factory (this is a test program specific for DWF)
  QpropDWF* QuarkPropDW = static_cast<QpropDWF*>(QP_DomainWallFact.getQuarkProp(conf_));
  // the prevoius static_cast is absolutely safe since we know exaclty what class we are creating
  
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

  ////////// test of Action_Nf2_ratio::DdagD2_inv  
  double mq1 = 0.05;
  double mq2 = 0.10;
  
  Dirac_optimalDomainWall Ddwf1(b,c,M0,mq1,omega,&(conf_.U));
  Dirac_optimalDomainWall Ddwf2(b,c,M0,mq2,omega,&(conf_.U));

  Fopr_DdagD DdagD2(&Ddwf2);

  Solver_CG SolvR2(10e-12,1000,&DdagD2);
  Field  D1phi =Ddwf1.mult_dag(phi);
  double phisum= D1phi.norm();
  double phinorm= Communicator::instance()->reduce_sum(phisum);

  CCIO::cout<<"phi_.norm="<<sqrt(phinorm)<<std::endl;
  CCIO::cout<<"typeid(D1_)="<<typeid(Ddwf1).name()<<std::endl;

  SolverOutput monitor;
  double diff;
  Field sol(Ddwf1.fsize());
  monitor = SolvR2.solve(sol,D1phi);
  monitor.print();
  ///////////////

  
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
  mult5d_gamma5_test(*DiracODWF,phi,100);
  CCIO::cout<< "Time for 100 calls : " 
      << ( ( clock() - start ) / (double)CLOCKS_PER_SEC ) <<'\n';
  
  // quark propagator
  double stop_cond = 1.0e-24;
  int    Niter= 10000;

  // It follows a standard construction (factories will use a similar one)
  //Dirac_optimalDomainWall Ddwf_PV(Ddwf_5d, PauliVillars);
  Dirac_optimalDomainWall Ddwf_PV(*DiracODWF, PauliVillars);
  // Solver* SolvDWF = new Solver_CG(stop_cond,Niter,new Fopr_DdagD(&Ddwf_5d));
  Solver* SolvDWF = new Solver_CG(stop_cond,Niter , new Fopr_DdagD(DiracODWF));
  Solver* SolvPV  = new Solver_CG(stop_cond,Niter , new Fopr_DdagD(&Ddwf_PV ));
  Dirac_optimalDomainWall_4D DiracDWF_4d(*DiracODWF, SolvDWF, SolvPV);
  QpropDWF QuarkPropagator(DiracDWF_4d);
  //////////////////////////////////// 
 
  CCIO::cout << ".::: Test Dirac_optimalDomainWall meson correlator" 
	     <<std::endl;
  // Here uses the default constructor with default solver
  //QpropDWF QuarkPropagator(Ddwf_5d,stop_cond,Niter);
  
  vector<int> spos(4,0); 
  Source_local<Format::Format_F> src(spos,CommonPrms::instance()->Nvol());
  
  prop_t sq;
  //QuarkPropagator.calc(sq,src);
  QuarkPropDW->calc(sq,src);

  
  GammaMatrices::Unit Gamma;
  MesonCorrelator meson(Pion);
  vector<double> mcorr = meson.calculate<Format::Format_F>(sq,sq);
  vector<double>::const_iterator it=mcorr.begin();
  int t=0;
  while(it!=mcorr.end()) CommunicatorItems::pprintf ("%d %.8e\n",t++, *it++);
  

  return 0;
}
