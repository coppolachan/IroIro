/*!
 * This still needs the routine to save the solution vectors
 * phi(x,t)^{(i)}_[r]
 *   (i is dilution index and r is the random number label
 *    as described in the notes)
 * and the routines to read in the scalarfield_3D
 *
 * This routine will only construct the solutions for r=0
 * (as in a single random source) 
 * For 4 random sources, one must repeat this calculation with 
 * a different random number seed
 * 
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
using namespace EvenOddUtils;

typedef FermionField1sp ScalarField;

int Test_LapH_Solver::run()
{
  CCIO::cout<<"Test_LapH_Solver::run() called\n";

  // Lattice size parameters
  int Nvol = CommonPrms::instance()->Nvol();

  int Nt = CommonPrms::instance()->Nt();
  int Nx = CommonPrms::instance()->Nx();
  int Ny = CommonPrms::instance()->Ny();
  int Nz = CommonPrms::instance()->Nz();
  int Nvol3D = Nx*Ny*Nz;
 
  Staples stpl;
  double plq = stpl.plaquette(Gfield_);
  CCIO::cout<<" Plaquette ="<< plq <<std::endl;
  
  /************************************************************************************/
  //
  // For 5-D Inversion (from test_DWF.cpp)
  //
  Field solution;
  Field* u = &(Gfield_.data);
  double M0 = -1.0;
  // creation of Dirac_Wilson operators 
  Dirac_Wilson Dw_eo(M0,u,Dop::EOtag());
  Dirac_Wilson Dw_oe(M0,u,Dop::OEtag());
  // Definition of the 5D DWF object
  int N5=4;
  double b=2.0;
  double c=1.0;
  //###########################################################
  double mq =0.007; // one should read this from the xml as well 
  //###########################################################

  int Nev = 120; // should read these in from xml ############################################
  int Nevdil = 6; //interlace-6
  int Nspindil = ND_; // full spin dilution
  int Ndil = Nt * Nevdil * Nspindil; // 64 * 6   * 4 =  1536;// max dilution index + 1 (but we'll only use t=0)


  ffmt_t fmt5d(CommonPrms::instance()->Nvol(),N5);
  // For the TanH approximation
  std::vector<double> omega(N5,1.0);
 
  // Solver parameters
  double prec = 1.0e-16;
  int max_iter = 600;
  
  // Constructs the 5D objects for the operator and the Pauli-Villars
  Dirac_DomainWall_EvenOdd Ddwf_eo(b,c,M0, mq,omega,&Dw_eo,&Dw_oe);
  Dirac_DomainWall_EvenOdd Ddpv_eo(b,c,M0,1.0,omega,&Dw_eo,&Dw_oe);
  Fopr_DdagD DdagDdwf_eo(&Ddwf_eo);
  Fopr_DdagD DdagDdpv_eo(&Ddpv_eo);
  Solver_CG slv_dwf_eo(prec,max_iter,&DdagDdwf_eo);
  Solver_CG slv_dpv_eo(prec,max_iter,&DdagDdpv_eo); 
  Inverter_WilsonLike invDdwf(&Ddwf_eo,&slv_dwf_eo);
  Inverter_WilsonLike invDdpv(&Ddpv_eo,&slv_dpv_eo);
  Dirac_DomainWall_4D_eoSolv D4eo(N5,mq,&invDdwf,&invDdpv);
  
  /*************************************************************************/
  //
  // Generating The Source Vector 
  //
  //std::vector< ScalarField > Sfield;
  ScalarField Sfield(Nvol3D); // Create a 3d object
  // Creating the 4d fermion field (color and spin dof)
  
  FermionField Ffield;
  FermionField srcfield;
  
  //
  // Eigenvectors are stored as 3-D vectors for each timeslice
  // i.e. "eigen_t0_5830" contains 120, 3D-eigenvectors on timeslice t=0 for config 5830
  //     
  // So at this stage, we need to read in 120 3-D eigenvectors from the file "eigen_t0_####"
  // into a 3-D field, Sfield_3D, I think. It's the same if t0=0
  //
  // For now, we will just need the t=t0=0 eigenvectors ...
  //

  //int num_rand = Nt * Nev * Nspindil;// 64 * 120 * 4 = 30720;// number of random numbers to generate
  int num_rand =  Nev * Nspindil;// 64 * 120 * 4 = 7200;// number of random numbers to generate

  // #######################################################################################
  // this one needs some redesign to deal with 3d vectors
  // now loading 8 (64/8 <- number of nodes in the T direction) 
  //       The filename should come from the xml file 

  // #######################################################################################
  // 
  // Creating the RNG from the XML file
  RNG_Env::RNG = RNG_Env::createRNGfactory(LapH_node);
  RandNum* rand_ = RNG_Env::RNG->getRandomNumGenerator();
  //
  // Example:
  //  Connected diagrams  Tf, Gf, EVi6
  //  Nt=64  Tf   full-time  (32^3 x 64 lattice)
  //  ND=4   Gf   full-spin
  //  Nev=6  EVi6 interlace-6  (120 eigenvectors of 3-D laplacian)
  //      Need: 120*4 = 480  Z4 random numbers at each t0
  //

  std::valarray<double> white_noise(0.0, num_rand);
  std::valarray<double> rho(0.0, num_rand*2); 
  rand_->get(white_noise);    
  // Make Z4 noise from uniform noise
  for (int idx = 0; idx <white_noise.size(); idx++)
  {
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
  // *****************************************************************************
  //
  // Finding the solutions for t0=0 source only for now
  // (connected diagrams only)
  //
  // Disconnected diagrams will be dealt with later
  //
  int time_slice = 0;
  //
  // rho has ev, spin and t indices
  //    (here, t is essentially absent from the indices because timeslice t=0)
  //
  int rho_idx;  
  for(int ev_start=0; ev_start < Nevdil; ++ev_start)
    {
      for(int eigvec_number=ev_start; eigvec_number < Nev; eigvec_number+=Nevdil)
	{     
	  CCIO::ReadFromDisk3D< Format::Format_S >(Sfield.data, "eigen_t0_5830", eigvec_number, 0, "Binary", false);
	  srcfield = 0;
	  for(int spin=0; spin < ND_; ++spin)
	    {
	    
	      
	      rho_idx = (ND_* eigvec_number + spin); // select random number
	      
	      ScalarField Sfield_3D( Nvol3D );
	      //copy the #eigenvec_number vector into Sfield_3D
	      for (int i = 0; i  < Sfield_3D.size(); i++) 
		{
		  // this assumes t0 = 0                              check Nvol
		  Sfield_3D.data.set(i, Sfield.data[i]);
		}
	      // Dilute the spin ND_ = 4  NC_ = 3 Ffield.format.Nin() = 2*NC_*ND_ number of real dof
	      int offset = time_slice*(Ffield.format.Nin()*Nvol3D)+2*NC_*spin;
	      for (int i = 0; i  < Sfield_3D.size(); i+=2*NC_) 
		{
		  for (int s = 0; s < 2*NC_ ; s++ ) 
		    Ffield.data.set(offset + i*ND_ + s, Sfield_3D.data[i+s]);
		}
	      // Multiply by the random number with index rho_idx
	      for (int  i=0; i  < Ffield.size(); i+=2) {
		double real = rho[2*rho_idx]*Ffield[i] - rho[2*rho_idx+1]*Ffield[i+1];
		double imag = rho[2*rho_idx]*Ffield[i+1] + rho[2*rho_idx+1]*Ffield[i];
		Ffield.data.set(i, real);
		Ffield.data.set(i, imag);
	      }
	      // add to the source field
	      srcfield += Ffield;    
	    }//spin
	} // evigvec_number
	  /*
      for (int i = 0; i  < srcfield.size(); i+=2) {
	CCIO::cout << "["<<i<<"] F = "<< srcfield[i] << " , " << srcfield[i+1] << "\n";
      }
      */

      // Save the 4-D source field for now: rho_{i}_s{spin}_t{t0}_noise_0(x,t) for noise 0
      char filename[128]; // Need input from xml ##############################################
      
      
      // Invert to find the solution
      sprintf(filename, "smeared_rho_0_NvI6_t0_i%d", ev_start);
      CCIO::SaveOnDisk<FermionField>(srcfield.data, filename, true);
      Field solution = D4eo.mult_inv(srcfield.data);  
	
      
      // Save the 4-D fermion solution: phi_{i}_t{t0}_noise_0(x,t) 
      sprintf(filename, "phi_0_NvI6_t0_i%d", ev_start);
      CCIO::SaveOnDisk<FermionField>(solution, filename, true);
      // One can also project the solution onto the eigenvectors at this stage and
      // save only the coefficients if we know we will be smearing at the sink.   
      /*
      for(int site=0; site<Nvol; ++site){
	double Op = 0;
	if (SiteIndex::g_t(site) == time_slice){
	  Op += srcfield.data[srcfield.format.slice_islice(site)]*(Dw_eo.gamma5(solution));
	}
      }
      */
      double Op = srcfield.data*(Dw_eo.gamma5(solution));

    } // ev_start

  return 0;
}

