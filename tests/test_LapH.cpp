/*!
 * @file test_LapH.cpp
 * @brief Tests for the LapH solver
 */
#include <stdio.h>

#include "test_LapH.hpp"
#include "Dirac_ops/dirac_DomainWall.hpp"
#include "Dirac_ops/dirac_DomainWall_EvenOdd.hpp"
#include "Dirac_ops/dirac_DomainWall_4D_eoSolv.hpp"
#include "Dirac_ops/dirac_DomainWall_4D_fullSolv.hpp"
#include "Solver/solver_CG.hpp"
#include "Solver/solver_BiCGStab.hpp"
#include "Tools/RandomNumGen/randNum_MT19937.h"
#include "Measurements/GaugeM/staples.hpp"

using namespace std;
using namespace Format;
using namespace DomainWallFermions;
using namespace EvenOddUtils;

typedef FermionField1sp ScalarField;

int Test_LapH_Solver::run(){
  CCIO::cout<<"Test_LapH_Solver::run() called\n";
  RNG_Env::RNG = RNG_Env::createRNGfactory(LapH_node);

  int Nvol = CommonPrms::instance()->Nvol();
  Staples stpl;
  double plq = stpl.plaquette(Gfield_);
  CCIO::cout<<" plaq="<<plq<<std::endl;

  ///// Generating source vector 
  //std::vector< ScalarField > Sfield;
  ScalarField Sfield;

  //CCIO::ReadFromDisk< Format::Format_S >((Sfield[0]).data, "eigen_t0_5830", 1);
  CCIO::ReadFromDisk< Format::Format_S >(Sfield.data, "eigen_t0_5830", 0, "Binary", false);


  RandNum* rand_ = RNG_Env::RNG->getRandomNumGenerator();

  FermionField Ffield;
  int Nx = CommonPrms::instance()->Nx();
  int Ny = CommonPrms::instance()->Ny();
  int Nz = CommonPrms::instance()->Nz();

  /////////////////////////////////
  int spin = 0;
  int time_slice = 0;
  /////////////////////////////////

  //////////////////////////////////////
  int num_rand = 4;
  //////////////////////////////////////

  std::valarray<double> white_noise(0.0, num_rand);
  std::valarray<double> rho(0.0, num_rand*2);

    
  rand_->get(white_noise);
    
  for (int idx = 0; idx <white_noise.size(); idx++){
    if ( white_noise[idx] < 0.25){
    rho[2*idx]   = 1.0;
    rho[2*idx+1] = 0.0;        
    continue;
    }
    if ( white_noise[idx] < 0.5){
    rho[2*idx]   = -1.0;
    rho[2*idx+1] = 0.0;        
    continue;
    }
    if ( white_noise[idx] < 0.75){
    rho[2*idx]   = 0.0;
    rho[2*idx+1] = 1.0;        
    continue;
    }
    rho[2*idx]   = 0.0;
    rho[2*idx+1] = -1.0;        
  }



  ////////////////////////
  // Spin dilution
  int eigvec_number = 0;
  
  ScalarField Sfield_3D(Nvol/CommonPrms::instance()->Nt());
  //copy the first vector into Sfield_3D
  
  for (int i = 0; i  < Sfield_3D.size(); i++) {
    Sfield_3D.data.set(i, Sfield.data[i+Sfield.format.Nin()*Nvol*eigvec_number]);
  }


  int offset = time_slice*(Ffield.format.Nin()*Nx*Ny*Nz)+2*NC_*spin;
  for (int i = 0; i  < Sfield_3D.size(); i+=2*NC_) {
    Ffield.data.set(offset + i*ND_   , Sfield[i]);
    Ffield.data.set(offset + i*ND_+1 , Sfield[i+1]);
    Ffield.data.set(offset + i*ND_+2 , Sfield[i+2]);
    Ffield.data.set(offset + i*ND_+3 , Sfield[i+3]);
    Ffield.data.set(offset + i*ND_+4 , Sfield[i+4]);
    Ffield.data.set(offset + i*ND_+5 , Sfield[i+5]);
  
  }

  for (int i = 0; i  < Ffield.size(); i+=2) {
    CCIO::cout << "["<<i<<"] F = "<< Ffield[i] << " , " << Ffield[i+1] << "\n";
  }

  int rho_idx = 0;
  for (int  i = 0; i  < Ffield.size(); i+=2) {
    double real = rho[2*rho_idx]*Ffield[i] - rho[2*rho_idx+1]*Ffield[i+1];
    double imag = rho[2*rho_idx]*Ffield[i+1] + rho[2*rho_idx+1]*Ffield[i];
    Ffield.data.set(i, real);
    Ffield.data.set(i, imag);
  }




  /************************************************************************************/
  Field* u = &(Gfield_.data);
  double M0 = -1.0;

  // creation of Dirac_Wilson operators 
  Dirac_Wilson Dw_eo(M0,u,Dop::EOtag());
  Dirac_Wilson Dw_oe(M0,u,Dop::OEtag());
 
  // Definition of the 5D DWF object
  int N5=4;
  double b=2.0;
  double c=0.0;
  double mq =0.05;
  ffmt_t fmt5d(CommonPrms::instance()->Nvol(),N5);

  // For the TanH approximation
  std::vector<double> omega(N5,1.0);

  // Solver parameters
  double prec = 1.0e-16;
  int max_iter = 600;

  ///////
  // Constructs the 5D objects for the operator and the Pauli-Villars
  Dirac_optimalDomainWall_EvenOdd Ddwf_eo(b,c,M0, mq,omega,&Dw_eo,&Dw_oe);
  Dirac_optimalDomainWall_EvenOdd Ddpv_eo(b,c,M0,1.0,omega,&Dw_eo,&Dw_oe);

  Fopr_DdagD DdagDdwf_eo(&Ddwf_eo);
  Fopr_DdagD DdagDdpv_eo(&Ddpv_eo);

  Solver_CG slv_dwf_eo(prec,max_iter,&DdagDdwf_eo);
  Solver_CG slv_dpv_eo(prec,max_iter,&DdagDdpv_eo);

  Inverter_WilsonLike invDdwf(&Ddwf_eo,&slv_dwf_eo);
  Inverter_WilsonLike invDdpv(&Ddpv_eo,&slv_dpv_eo);

  Dirac_optimalDomainWall_4D_eoSolv D4eo(N5,mq,&invDdwf,&invDdpv);

  /*--------mult test---------*/

 // e/o 5D DWF
  Field Sdweo(Ddwf_eo.fsize()); 
  if(Communicator::instance()->primaryNode()) Sdweo.set(0,1.0);
  if(Communicator::instance()->primaryNode()) Sdweo.set(5,1.0);
  if(Communicator::instance()->primaryNode()) Sdweo.set(10,1.0);
  Field Wdweo = Ddwf_eo.mult(Sdweo);
  double Ndweo =Wdweo.norm();
  CCIO::cout<<"Ndweo="<<Ndweo<<"\n";

  double nwfl=0.0; 
  double nweo=0.0; 

  CCIO::cout<<"test of 4D e/o solver \n";
  Field S4eo(D4eo.fsize()); 
  if(Communicator::instance()->primaryNode()) S4eo.set(0,1.0);
  if(Communicator::instance()->primaryNode()) S4eo.set(5,1.0);
  if(Communicator::instance()->primaryNode()) S4eo.set(10,1.0);
  Field weo = D4eo.mult(S4eo); 
  nweo = weo.norm();

  CCIO::cout<<"nwfl="<<nwfl<<" nweo="<<nweo<<"\n";

  //////////////////////////////////////////////////////
  // Test the solver
  Field solution = D4eo.mult_inv(Ffield.data);


  // Solution should be stored on the disk

  return 0;
}

