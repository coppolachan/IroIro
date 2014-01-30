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
#include "Dirac_ops/dirac_Operator_FactoryCreator.hpp"
#include "Measurements/GaugeM/staples.hpp"

#include "include/timings.hpp"
#include "include/messages_macros.hpp"

using namespace std;
using namespace Format;
using namespace EvenOddUtils;

typedef FermionField1sp ScalarField;

int Test_LapH_Solver::run()
{
  CCIO::cout<<"Test_LapH_Solver::run() called\n";

  XML::node LapH_node = input_.node;
  InputConfig config = input_.getConfig();

  XML::descend(LapH_node,"LapH",MANDATORY);


  // Lattice size parameters
  int Nvol = CommonPrms::instance()->Nvol();
  int Nt = CommonPrms::instance()->Nt();
  int Nx = CommonPrms::instance()->Nx();
  int Ny = CommonPrms::instance()->Ny();
  int Nz = CommonPrms::instance()->Nz();
  int Nvol3D = Nx*Ny*Nz;
  Field solution;
  Field* u = &((input_.gconf)->data);

  Staples stpl;
  double plq = stpl.plaquette(*input_.gconf);
  CCIO::cout<<" Plaquette ="<< plq <<std::endl;
  
  /************************************************************************************/
  //
  // For 5-D Inversion (from test_DWF.cpp)
  //

  XML::descend(LapH_node,"KernelDWF_4d");
  auto_ptr<DiracDWF4dFactory> 
    Wilson_Kernel_4d_factory(Diracs::createDiracDWF4dFactory(LapH_node));
  auto_ptr<Dirac_DomainWall_4D> Wilson_Kernel_4d(Wilson_Kernel_4d_factory->getDirac(config));


  /*
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
  */

  int Nev = 120; // should read these in from xml ############################################
  int Nevdil = 6; //interlace-6
  int Nspindil = ND_; // full spin dilution
  int Ndil = Nt * Nevdil * Nspindil; // 64 * 6   * 4 =  1536;// max dilution index + 1 (but we'll only use t=0)

  /*
  ffmt_t fmt5d(CommonPrms::instance()->Nvol(),N5);
  // For the TanH approximation
  std::vector<double> omega(N5,1.0);
 
  // Solver parameters
  double prec = 1.0e-16;
  int max_iter = 5000;
  
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
  */
  /*************************************************************************/
  //
  // Generating The Source Vector 
  //
  //std::vector< ScalarField > Sfield;
  ScalarField Sfield(Nvol3D); // Create a 3d object
  // Creating the 4d fermion field (color and spin dof)
  
  FermionField Ffield;
  std::vector<FermionField> srcfield(ND_);
  
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


  // 
  // Creating the RNG from the XML file

  //RNG_Env::RNG = RNG_Env::createRNGfactory(LapH_node);
  //RandNum* rand_ = RNG_Env::RNG->getRandomNumGenerator();

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
  input_.rng->get(white_noise);    
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
  
  long double timer_source, timer_solver, timer_save;

  for(int ev_start=0; ev_start < Nevdil; ++ev_start)
    {
      FINE_TIMING_START(timer_source);


      for(int spin=0; spin < ND_; ++spin)
	srcfield[spin] = 0;

      for(int eigvec_number=ev_start; eigvec_number < Nev; eigvec_number+=Nevdil)
	{     
	  CCIO::ReadFromDisk3D< Format::Format_S >(Sfield.data, "eigen_t0_5830", eigvec_number, 0, "Binary", false);

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
	      srcfield[spin] += Ffield;    
	    }//spin
	} // evigvec_number

      FINE_TIMING_END(timer_source);
      CCIO::cout << "[Timing] Source(s) creation: "<< timer_source << "\n";
	  /*
      for (int i = 0; i  < srcfield.size(); i+=2) {
	CCIO::cout << "["<<i<<"] F = "<< srcfield[i] << " , " << srcfield[i+1] << "\n";
      }
      */

      // Save the 4-D source field for now: rho_{i}_s{spin}_t{t0}_noise_0(x,t) for noise 0
      char filename[128]; // Need input from xml ##############################################
      
      
      // Invert to find the solution
      for (int spin = 0; spin < ND_; spin++){
	CCIO::cout << " Ev_start: "<< ev_start << "  Spin : "<<spin << "\n";
	
	sprintf(filename, "smeared_rho_0_NvI6_t0_i%d", ev_start*ND_+spin);
	FINE_TIMING_START(timer_save);
	CCIO::SaveOnDisk<FermionField>(srcfield[spin].data, filename);
	FINE_TIMING_END(timer_save);
	CCIO::cout << "[Timing] Source(s) save: "<< timer_save << "\n";
	FINE_TIMING_START(timer_solver);
	Field solution = Wilson_Kernel_4d->mult_inv(srcfield[spin].data);  
  	FINE_TIMING_END(timer_solver);
      	CCIO::cout << "[Timing] Solver: "<< timer_solver << "\n";

	// Save the 4-D fermion solution: phi_{i}_t{t0}_noise_0(x,t) 
  	sprintf(filename, "phi_0_NvI6_t0_i%d", ev_start*ND_+spin);
	FINE_TIMING_START(timer_save);
	CCIO::SaveOnDisk<FermionField>(solution, filename);
	FINE_TIMING_END(timer_save);
	CCIO::cout << "[Timing] Solution(s) save: "<< timer_save << "\n";
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
      }
      

    } // ev_start

  return 0;
}

