/*!
 * @file dirac_BFM_wrapper_mixed_prec.cpp
 * @brief Defines the wrapper classs for P. Boyle Bagel/BFM libs - Mixed precision libraries
 * Time-stamp: <2014-04-25 13:09:57 neo>
 */

#include "dirac_BFM_wrapper.hpp"

#include "include/timings.hpp"
#include "include/messages_macros.hpp"
#include "include/errors.hpp"
#include <stdlib.h>     /* atoi */
#include <stdio.h>
#include <omp.h>


int Dirac_BFM_Wrapper::CGNE_mixed_prec(Fermion_t sol, Fermion_t src, int max_outer, double shift){
  ////////////////////
  int cb= Even;
  int iters = 0;
  int converged = 0;

  Fermion_t mtmp = linop.threadedAllocFermion();
  Fermion_t mmtmp = linop.threadedAllocFermion();
  
  Fermion_t tmp_d[2];
  Fermion_t tmp_s[2];
  Fermion_t sol_s[2];
  
  // Allocate internal fermions
  for(int cb=0;cb<2;cb++){
    tmp_d[cb]    = linop.threadedAllocFermion();
    tmp_s[cb]    = linop_single.threadedAllocFermion();
    sol_s[cb]    = linop_single.threadedAllocFermion();
  }
  
  int me = linop.thread_barrier();
  double defect_norm = 1.0;
  double true_residual;
  double src_norm = sqrt(linop.norm(src));
  
  for (int outer = 0; outer < max_outer; outer++){
    int inner_iters;
    if ( !me ) linop.comm_init();  // Start double comms
    linop.thread_barrier();
    
    linop.Mprec(sol, mtmp, tmp_d[0], 0);
    linop.Mprec(mtmp, mmtmp, tmp_d[0], 1);
    if (shift !=0.0) {
      linop.axpy(mtmp,sol,mmtmp,shift);
      linop.axpy(tmp_d[0], mtmp, src, -1.0);
    } else {
    linop.axpy(tmp_d[0], mmtmp, src, -1.0); // calculate defect (tmp_d) double prec
    }
    defect_norm = sqrt(linop.norm(tmp_d[0]));
    
    true_residual = defect_norm/src_norm;// since initial guess is zero
    
    if ( linop.isBoss() && (!me) ) printf("solve_CGNE_mixed_prec[%d] - defect norm  : %le\n",outer,defect_norm);
    if ( linop.isBoss() && (!me) ) printf("solve_CGNE_mixed_prec[%d] - true residual: %le\n",outer,true_residual);
    
    if ( !me ) linop.comm_end(); // End single comms
    linop.thread_barrier();
    
    
    // Check convergence
    if ( true_residual < linop.residual ) {
      if ( linop.isBoss() && !me ) 
	printf("solve_CGNE_mixed_prec[%d] - Iterations = %d  Residual = %le\n",
	       outer,iters,true_residual);

      outer = max_outer;
      converged=1;
    }
    
    if ( !converged ) { 
      
      //scale defect (to make efficient use of all lower precision digits)
      linop.scale(tmp_d[0], 1.0/defect_norm);
      
      ////////////////////////////////////////////////////
      // Single precision inner CG
      ////////////////////////////////////////////////////
      for(int cb=0;cb<2;cb++){
	//convert source from double to single
	linop.precisionChange(tmp_d[cb],tmp_s[cb],DoubleToSingle,cb);
	// Initial guess is set to zero. 
	linop_single.set_zero(sol_s[cb]);
      }
      linop.thread_barrier();
      
      
      double target_residual = linop.residual/true_residual / 10;
      if ( target_residual < 1.0e-6 ) target_residual = 1.0e-6;
      if ( linop_single.isBoss() && (!me) )
	printf("solve_CGNE_mixed_prec[%d]: setting target residual to %le\n",outer,target_residual);
      linop_single.residual = target_residual;
      
      if ( !me ) linop_single.comm_init();  // Start double comms
      linop_single.thread_barrier();
      
      // Iterate inner solver until ||b - Ax|| < linop_single.residual
      if (shift ==0){
      inner_iters = linop_single.CGNE_prec(sol_s[cb],tmp_s[cb]);
      } else {
	inner_iters = linop_single.CGNE_single_shift(sol_s[cb],tmp_s[cb], shift);
      }
      if( me == 0 ) {
	iters += inner_iters;
      }	  
      
      if ( !me ) linop_single.comm_end(); // End single comms
      linop_single.thread_barrier();
      
      ////////////////////////////////////////////////////
      // Convert to double and update solution
      ////////////////////////////////////////////////////
      // chi_h_{k+1}    = chi_h_k + defect_norm * (sol_s->double)
      for(int cb=0;cb<2;cb++){
	linop.precisionChange(sol_s[cb],tmp_d[cb],SingleToDouble,cb);// approx solution
	linop.axpy(sol,tmp_d[cb], sol,defect_norm);//new solution
      }
      
    }
    
  }

  
  //free allocated resources
  for(int cb=0;cb<2;cb++){
    linop.threadedFreeFermion(tmp_d[cb]);
    linop_single.threadedFreeFermion(tmp_s[cb]);
    linop_single.threadedFreeFermion(sol_s[cb]);
  }
  linop.threadedFreeFermion(mtmp);
  linop.threadedFreeFermion(mmtmp);
  linop.thread_barrier();
  ////////////////////////
  return converged;
}


int Dirac_BFM_Wrapper::CGNE_mixed_prec_residual_guided(Fermion_t sol, Fermion_t src, int max_outer){
  ////////////////////
  int cb= Even;
  int iters = 0;
  int converged = 0;

  Fermion_t mtmp = linop.threadedAllocFermion();
  Fermion_t mmtmp = linop.threadedAllocFermion();
  
  Fermion_t tmp_d[2];
  Fermion_t tmp_s[2];
  Fermion_t sol_s[2];
  Fermion_t src_s;
  Fermion_t p_s;

  // Allocate internal fermions
  for(int cb=0;cb<2;cb++){
    tmp_d[cb]    = linop.threadedAllocFermion();
    tmp_s[cb]    = linop_single.threadedAllocFermion();
    sol_s[cb]    = linop_single.threadedAllocFermion();
    src_s        = linop_single.threadedAllocFermion();
    p_s          = linop_single.threadedAllocFermion();
  }
  linop.precisionChange(src,src_s,DoubleToSingle,cb);


  int me = linop.thread_barrier();
  double defect_norm = 1.0;
  double defect_old;
  double true_residual;
  double src_norm = sqrt(linop.norm(src));
  
  double alpha = 0.0;
  double beta  = 0.0;
  double rho;

  ///////// initialization of the outer loop
  if ( !me ) linop.comm_init();  // Start double comms
  linop.thread_barrier();
  
  linop.Mprec(sol, mtmp, tmp_d[0], 0);
  linop.Mprec(mtmp, mmtmp, tmp_d[0], 1);
  linop.axpy(tmp_d[0], mmtmp, src, -1.0); // calculate defect (tmp_d) double prec
  

  defect_norm = sqrt(linop.norm(tmp_d[0]));
  
  true_residual = defect_norm/src_norm;// since initial guess is zero
  
  if ( linop.isBoss() && (!me) ) printf("solve_CGNE_mixed_prec init - defect norm  : %le\n",defect_norm);
  if ( linop.isBoss() && (!me) ) printf("solve_CGNE_mixed_prec init - true residual: %le\n",true_residual);
  
  if ( !me ) linop.comm_end(); // End single comms
  linop.thread_barrier();
 
  //scale defect (to make efficient use of all lower precision digits)
  linop.scale(tmp_d[cb], 1.0/defect_norm);
  //linop_single.scale(src_s, 1.0/defect_norm);
  /////////////////////////////////////
  linop.precisionChange(tmp_d[cb],p_s,DoubleToSingle,cb);
  ////////////////////////////////


  for (int outer = 0; outer < max_outer; outer++){
    int inner_iters;

    if ( !converged ) { 
      
      ////////////////////////////////////////////////////
      // Single precision inner CG
      ////////////////////////////////////////////////////
      for(int cb=0;cb<2;cb++){
	//convert source from double to single
	linop.precisionChange(tmp_d[cb],tmp_s[cb],DoubleToSingle,cb);
	// Initial guess is set to zero. 
	linop_single.set_zero(sol_s[cb]);
      }
      linop.thread_barrier();
      
      if (outer !=0 ){
	double proj = linop_single.inner_real(p_s, tmp_s[cb]);
	linop_single.axpy(p_s, tmp_s[cb], p_s, -proj);
	linop_single.axpy(p_s, p_s, tmp_s[cb], beta);
      }
 
      //linop_single.max_iter = 50;

      if ( !me ) linop_single.comm_init();  // Start double comms
      linop_single.thread_barrier();
      
      // Iterate inner solver until ||b - Ax|| < linop_single.residual
      inner_iters = linop_single.CGNE_prec_inner(sol_s[cb],src_s, tmp_s[cb], p_s, 1e-8);
 

      if( me == 0 ) {
	iters += inner_iters;
      }	  
      
      if ( !me ) linop_single.comm_end(); // End single comms
      linop_single.thread_barrier();
      
      ////////////////////////////////////////////////////
      // Convert to double and update solution
      ////////////////////////////////////////////////////
      // chi_h_{k+1}    = chi_h_k + defect_norm * (sol_s->double)
      rho = linop_single.norm(tmp_s[cb]);
      for(int cb=0;cb<2;cb++){
	linop.precisionChange(sol_s[cb],tmp_d[cb],SingleToDouble,cb);// approx solution
	linop.axpy(sol,tmp_d[cb], sol,defect_norm);//new solution
      }
      
      
    }
    
    if ( !me ) linop.comm_init();  // Start double comms
    linop.thread_barrier();
    
    linop.Mprec(sol, mtmp, tmp_d[0], 0);
    linop.Mprec(mtmp, mmtmp, tmp_d[0], 1);
    linop.axpy(tmp_d[0], mmtmp, src, -1.0); // calculate defect (tmp_d) double prec
  
    defect_old = defect_norm;
    defect_norm = sqrt(linop.norm(tmp_d[0]));
    beta = defect_norm/(defect_old*rho);


    true_residual = defect_norm/src_norm;// since initial guess is zero
    
    if ( linop.isBoss() && (!me) ) printf("solve_CGNE_mixed_prec[%d] - defect norm  : %le\n",outer,defect_norm);
    if ( linop.isBoss() && (!me) ) printf("solve_CGNE_mixed_prec[%d] - true residual: %le\n",outer,true_residual);
    
    if ( !me ) linop.comm_end(); // End single comms
    linop.thread_barrier();
    
    //scale defect (to make efficient use of all lower precision digits)
    linop.scale(tmp_d[0], 1.0/defect_norm);
    //linop_single.scale(src_s, defect_old/defect_norm);
    linop_single.scale(p_s, 1.0/beta);   
    // Check convergence
    if ( true_residual < linop.residual ) {
      if ( linop.isBoss() && !me ) 
	printf("solve_CGNE_mixed_prec[%d] - Iterations = %d  Residual = %le\n",
	       outer,iters,true_residual);

      outer = max_outer;
      converged=1;
    }


  }

  
  //free allocated resources
  for(int cb=0;cb<2;cb++){
    linop.threadedFreeFermion(tmp_d[cb]);
    linop_single.threadedFreeFermion(tmp_s[cb]);
    linop_single.threadedFreeFermion(sol_s[cb]);
  }
  linop.threadedFreeFermion(mtmp);
  linop.threadedFreeFermion(mmtmp);
  linop.thread_barrier();
  ////////////////////////
  return converged;
}




void Dirac_BFM_Wrapper::solve_CGNE_mixed_prec(FermionField& solution, FermionField& source){
  // Reference papers
  // D. Göddeke,  Robert Strzodka, and Stefan Turek. 
  // Accelerating double precision FEM simulations with GPUs. 
  // In Proceedings of ASIM 2005 - 18th Symposium on Simulation Technique, Sep 2005
  //
  // Robert Strzodka and Dominik Göddeke. 
  // Pipelined mixed precision algorithms on FPGAs for fast and accurate PDE solvers from low precision components. 
  // In IEEE Symposium on Field-Programmable Custom Computing Machines (FCCM 2006), pages 259–268, April 2006.

  if(is_initialized){
    long double export_timing;
    long double solver_timing;

    // Load configurations 
    FINE_TIMING_START(export_timing); 
    BFM_interface.GaugeExport_to_BFM(u_);
    BFM_interface_single.GaugeExport_to_BFM(u_);
    FINE_TIMING_END(export_timing);
    _Message(TIMING_VERB_LEVEL, "[Timing] - Dirac_BFM_Wrapper::solve_CGNE_mixed prec"
	     << " - Gauge Export to BFM timing = "
	     << export_timing << std::endl); 

    int cb = Even;// by default
    
    // Load the fermion field to BFM
    // temporary source vector for the BFM data 
    FermionField BFMsource(Nvol_*BFMparams.Ls_);
    FermionField BFMsol(Nvol_*BFMparams.Ls_);

    int half_vec = source.size()/BFMparams.Ls_;
    for (int s =0 ; s< BFMparams.Ls_; s++){
      for (int i = 0; i < half_vec; i++){
	
	BFMsource.data.set(i+    2*s*half_vec, source.data[i+half_vec*s]);//just even part is enough
	BFMsource.data.set(i+(2*s+1)*half_vec, 0.0);//just even part is enough

	BFMsol.data.set(i+    2*s*half_vec, solution.data[i+half_vec*s]);//just even part is enough
	BFMsol.data.set(i+(2*s+1)*half_vec, 0.0);//just even part is enough
	

      }
    }
    LoadSource(BFMsource, cb);
    LoadGuess(BFMsol, cb);

    //////////////////////////////////
    int max_outer = linop.max_iter/100;
    int converged = 0;
    //////////////////////////// Execute the solver
    FINE_TIMING_START(solver_timing); 
    #pragma omp parallel
    {
    #pragma omp for 
      for (int t=0;t<threads;t++){
	converged = CGNE_mixed_prec(chi_h[cb], psi_h[cb], max_outer);// performs better
	//converged = CGNE_mixed_prec_residual_guided(chi_h[cb], psi_h[cb], max_outer);
      }
    }
    
    if (!converged) {
      ErrorString msg("BFM Mixed precision solver not converged");
      Errors::ConvergenceErr(msg);
    }
    
    FINE_TIMING_END(solver_timing);
    _Message(TIMING_VERB_LEVEL, "[Timing] - Dirac_BFM_Wrapper::solve_CGNE_mixed prec"
	     << " - Solver total timing = "
	     << solver_timing << std::endl); 

    ///////////////////////////////////////////////
    // Get the solution from BFM 
    GetSolution(BFMsol,cb);
    for (int s =0 ; s< BFMparams.Ls_; s++){
      for (int i = 0; i < half_vec; i++){
	solution.data.set(i+half_vec*s, BFMsol.data[i+2*s*half_vec]);//copy only the even part
      }
    }



  } else {
    CCIO::cout << "The operator was not initialized yet\n";
  }



}



void Dirac_BFM_Wrapper::solve_CGNE_multishift_mixed_precision(std::vector < FermionField > & solutions, 
							      FermionField& source,
							      vector_double shifts,
							      vector_double mresiduals){
  if(is_initialized){
    // auxiliary fields
    std::vector < FermionField > internal_sol(shifts.size());
    FermionField BFMsource(Nvol_*BFMparams.Ls_);
    Fermion_t ms_chi[shifts.size()];
    vector_double m_alpha(shifts.size());

    long double export_timing;
    long double solver_timing;
    int cb = Even;// by default
    int dont_sum = 0;
    int half_vec = source.size()/BFMparams.Ls_;

    // Load configurations 
    FINE_TIMING_START(export_timing); 
    BFM_interface.GaugeExport_to_BFM(u_);
    BFM_interface_single.GaugeExport_to_BFM(u_);
    FINE_TIMING_END(export_timing);
    _Message(TIMING_VERB_LEVEL, "[Timing] - Dirac_BFM_Wrapper::solve_CGNE_multishift_mixed precision"
	     << " - Gauge Export to BFM timing = "
	     << export_timing << std::endl); 

    // Allocate the solutions fields in BFM - double precision
    for (int shift = 0; shift < shifts.size(); shift++){
      solutions[shift].data.resize(source.size());      //half vector
      internal_sol[shift].data.resize(2*source.size()); //full vector
      ms_chi[shift] = linop.allocFermion();             //half vector
      m_alpha[shift] = 1.0; //fixed
    }

    // assumes even/odd indexing
    for (int s =0 ; s< BFMparams.Ls_; s++){
	memcpy((void*)&(BFMsource.data[2*s*half_vec]),
	       (void*)&(source.data[s*half_vec]), 
	       sizeof(double)*half_vec);
    }
    // Load the fermion field to BFM
    LoadSource(BFMsource, cb);


    //////////////////////////////////
    int max_outer = linop.max_iter/100;
    int converged = 0;
    //////////////////////////// Execute the solver
    FINE_TIMING_START(solver_timing); 
    #pragma omp parallel for
    for (int t=0;t<threads;t++){
      int iter = CGNE_prec_MdagM_multi_shift_mixed_prec(ms_chi,
							psi_h[cb],
							&shifts[0], 
							&m_alpha[0],
							shifts.size(), 
							&mresiduals[0],
							dont_sum,
							max_outer);
    }
    FINE_TIMING_END(solver_timing);
    _Message(TIMING_VERB_LEVEL, "[Timing] - Dirac_BFM_Wrapper::solve_CGNE_multishift_mixed_precision"
	     << " - Solver total timing = "
	     << solver_timing << std::endl); 
    
    // Get the solution from BFM and convert basis 
    GetMultishiftSolutions(internal_sol,ms_chi,cb);
    
    for (int shift = 0; shift < shifts.size(); shift++){
      for (int s =0 ; s< BFMparams.Ls_; s++){
	memcpy((void*)&(solutions[shift].data[half_vec*s]),
	       (void*)&(internal_sol[shift].data[2*s*half_vec]), 
	       sizeof(double)*half_vec);
      }
      linop.freeFermion(ms_chi[shift]);
    }
    
  } else {
    CCIO::cout << "The operator was not initialized yet\n";
  }
}

/*
int Dirac_BFM_Wrapper::CGNE_prec_MdagM_multi_shift_mixed_prec_test1(Fermion_t sol[], //solution
								    Fermion_t src,   //source
								    double    mass[],
								    double    alpha[],
								    int       nshift,
								    double mresidual[],
								    int single,
								    int max_outer){
  ////////////////////
  int cb= Even;
  int iters = 0;
  int converged = 0;
  int global_conv;
  //Save the mass 
  Fermion_t mtmp = linop.threadedAllocFermion();
  Fermion_t mmtmp = linop.threadedAllocFermion();
  
  Fermion_t tmp_d[nshift];
  Fermion_t tmp_s[2];
  Fermion_t sol_s[nshift];
  
  // Allocate internal fermions
  for (int shift = 0; shift < nshift; shift++){
    tmp_d[shift] = linop.threadedAllocFermion();
    sol_s[shift] = linop_single.threadedAllocFermion();
  }

  for(int cb=0;cb<2;cb++){
    tmp_s[cb]    = linop_single.threadedAllocFermion();
   
  }
  
  int me = linop.thread_barrier();
  double defect_norm[nshift];
  double true_residual[nshift];
  double target_res = mresidual[0];
  for (int shift = 0; shift < nshift; shift++){
    defect_norm[shift] = 1.0;
  }
  double src_norm = sqrt(linop.norm(src));
  
  for (int outer = 0; outer < max_outer; outer++){
    int inner_iters;
    if ( !me ) linop.comm_init();  // Start double comms
    linop.thread_barrier();
    
    global_conv = 1;
    // Calculate defect for every solution

    for (int shift = 0; shift < nshift; shift++){ 
      linop.Mprec(sol[shift], mtmp, tmp_d[shift], 0);
      linop.Mprec(mtmp, mmtmp, tmp_d[shift], 1);
      linop.axpy(mtmp,sol[shift],mmtmp,mass[shift]);
      linop.axpy(tmp_d[shift], mtmp, src, -1.0); // calculate defect (tmp_d[shift]) double prec
      
      defect_norm[shift] = sqrt(linop.norm(tmp_d[shift]));
      
      true_residual[shift] = defect_norm[shift]/src_norm;// since initial guess is zero
      
      if ( linop.isBoss() && (!me) ) printf("solve_CGNE_multishift_mixed_prec[%d] - defect norm[%d]    : %le\n",outer,shift,defect_norm[shift]);
      if ( linop.isBoss() && (!me) ) printf("solve_CGNE_multishift_mixed_prec[%d] - true residual[%d]  : %le\n",outer,shift,true_residual[shift]);
    }
    // Check convergence of lowest shift every solution
    if ( true_residual[0] < target_res ) {
      if ( linop.isBoss() && !me ) 
	printf("solve_CGNE_multishift_mixed_prec[%d] - Iterations = %d  Residual = %le\n",
	       outer,iters,true_residual[0]);
      
      outer = max_outer;
      converged=1;
      
    }
   
    
    
    if ( !me ) linop.comm_end(); // End single comms
    linop.thread_barrier();

    
    if ( !converged ) { 
      //scale defect (to make efficient use of all lower precision digits)
      for (int shift = 0; shift < nshift; shift++)
	linop.scale(tmp_d[shift], 1.0/defect_norm[shift]); //biggest defect
      
      ////////////////////////////////////////////////////
      // Single precision inner CG
      ////////////////////////////////////////////////////
      //convert source from double to single
      linop.precisionChange(tmp_d[0],tmp_s[cb],DoubleToSingle,cb);
      // Initial guess is set to zero. 
      for (int shift = 0; shift < nshift; shift++)
	linop_single.set_zero(sol_s[shift]);
      
      linop.thread_barrier();
      
      
      double target_residual = target_res/true_residual[0] / 10;
      if ( target_residual < 1.0e-6 ) target_residual = 1.0e-6;
      if ( linop_single.isBoss() && (!me) )
	printf("solve_CGNE_mixed_prec[%d]: setting target residual to %le\n",outer,target_residual);
      //linop_single.residual = target_residual;
      
      if ( !me ) linop_single.comm_init();  // Start double comms
      linop_single.thread_barrier();
      
      // Iterate inner solver until ||b - Ax|| < linop_single.residual
      // shifts must be reduced when converged
      for (int shift = 0; shift < nshift; shift++)
	mresidual[shift] = target_residual;
      inner_iters = linop_single.CGNE_prec_MdagM_multi_shift(sol_s,
							     tmp_s[cb],
							     &mass[0], 
							     &alpha[0],
							     nshift, 
							     &mresidual[0],
							     0);
      if( me == 0 ) {
	iters += inner_iters;
      }	  
      
      if ( !me ) linop_single.comm_end(); // End single comms
      linop_single.thread_barrier();
      
      ////////////////////////////////////////////////////
      // Convert to double and update solution
      ////////////////////////////////////////////////////
      // chi_h_{k+1}    = chi_h_k + defect_norm * (sol_s->double)
      for (int shift = 0; shift < nshift; shift++){
	linop.precisionChange(sol_s[shift],tmp_d[shift],SingleToDouble,cb);// approx solution
	linop.axpy(sol[shift],tmp_d[shift], sol[shift],defect_norm[shift]);//new solution
      }
    
    }
    
  }


  // Refine other solutions with mixed cg
  for (int s = 1; s < nshift; s++){
    CGNE_mixed_prec(sol[s], src, max_outer, mass[s]);
    
  }





  
  //free allocated resources
  for(int cb=0;cb<2;cb++){
    linop.threadedFreeFermion(tmp_d[cb]);
    linop_single.threadedFreeFermion(tmp_s[cb]);
    linop_single.threadedFreeFermion(sol_s[cb]);
  }
  linop.threadedFreeFermion(mtmp);
  linop.threadedFreeFermion(mmtmp);
  linop.thread_barrier();
  ////////////////////////
  return converged;
  
}
*/

int Dirac_BFM_Wrapper::CGNE_prec_MdagM_multi_shift_mixed_prec(Fermion_t sol[], //solution
							      Fermion_t src,   //source
							      double    mass[],
							      double    alpha[],
							      int       nshift,
							      double mresidual[],
							      int single,
							      int max_outer){
  ////////////////////
  int cb= Even;
  int iters = 0;
  int converged = 0;

  Fermion_t mtmp = linop.threadedAllocFermion();
  Fermion_t mmtmp = linop.threadedAllocFermion();
  
  Fermion_t tmp_d;
  Fermion_t tmp_s;
  Fermion_t sol_s[nshift];
  
  // Allocate internal fermions
  for (int shift = 0; shift < nshift; shift++){
    sol_s[shift] = linop_single.threadedAllocFermion();
  }

  tmp_d = linop.threadedAllocFermion();
  tmp_s = linop_single.threadedAllocFermion();
  
  int me = linop.thread_barrier();
  double defect_norm;
  double true_residual;
  double src_norm = sqrt(linop.norm(src));
  
  int inner_iters;
  if ( !me ) linop.comm_init();  // Start double comms
  linop.thread_barrier();
  
  // Calculate defect for smallest shift
  
  linop.Mprec(sol[0], mtmp, tmp_d, 0);
  linop.Mprec(mtmp, mmtmp, tmp_d, 1);
  linop.axpy(mtmp,sol[0],mmtmp,mass[0]);
  linop.axpy(tmp_d, mtmp, src, -1.0); // calculate defect (tmp_d[shift]) double prec
  
  defect_norm = sqrt(linop.norm(tmp_d));
  
  true_residual = defect_norm/src_norm;// since initial guess is zero
  
  if ( linop.isBoss() && (!me) ) printf("solve_CGNE_multishift_mixed_prec - defect norm    : %le\n",defect_norm);
  if ( linop.isBoss() && (!me) ) printf("solve_CGNE_multishift_mixed_prec - true residual  : %le\n",true_residual);
    
  // Check convergence of lowest shift every solution
  if ( true_residual < mresidual[0] ) {
    if ( linop.isBoss() && !me ) 
      printf("solve_CGNE_multishift_mixed_prec - Iterations = %d  Residual = %le\n",
	     iters,true_residual);
    converged=1;
  }
  
  if ( !me ) linop.comm_end(); // End single comms
  linop.thread_barrier();
  
  if ( !converged ) { 
    //scale defect (to make efficient use of all lower precision digits)
    linop.scale(tmp_d, 1.0/defect_norm); //biggest defect
    
    ////////////////////////////////////////////////////
    // Single precision inner CG
    ////////////////////////////////////////////////////
    //convert source from double to single
    linop.precisionChange(tmp_d,tmp_s,DoubleToSingle,cb);
    // Initial guesses is set to zero. 
    for (int shift = 0; shift < nshift; shift++)
      linop_single.set_zero(sol_s[shift]);
    
    linop.thread_barrier();
    
    
    double target_residual = mresidual[0]/true_residual / 10;
    if ( target_residual < 1.0e-6 ) target_residual = 1.0e-6;
    if ( linop_single.isBoss() && (!me) )
      printf("solve_CGNE_multishift_mixed_prec: setting target residual to %le\n",target_residual);
    //linop_single.residual = target_residual;
    
    if ( !me ) linop_single.comm_init();  // Start double comms
    linop_single.thread_barrier();
    
    // Iterate inner solver until ||b - Ax|| < linop_single.residual
    // shifts must be reduced when converged
    for (int shift = 0; shift < nshift; shift++)
      mresidual[shift] = target_residual;
    inner_iters = linop_single.CGNE_prec_MdagM_multi_shift(sol_s,
							   tmp_s,
							   &mass[0], 
							   &alpha[0],
							   nshift, 
							   &mresidual[0],
							   0);
    if( me == 0 ) {
      iters += inner_iters;
    }	  
    
    if ( !me ) linop_single.comm_end(); // End single comms
    linop_single.thread_barrier();
    
    ////////////////////////////////////////////////////
    // Convert to double and update solution
    ////////////////////////////////////////////////////
    // chi_h_{k+1}    = chi_h_k + defect_norm * (sol_s->double)
    for (int shift = 0; shift < nshift; shift++){
      linop.precisionChange(sol_s[shift],tmp_d,SingleToDouble,cb);// approx solution
      linop.axpy(sol[shift],tmp_d, sol[shift],defect_norm);//new solution
    }


  }
  
  // Refine solutions with mixed cg
  for (int s = 0; s < nshift; s++){
    CGNE_mixed_prec(sol[s], src, max_outer, mass[s]);
  }

  
  //free allocated resources
  linop.threadedFreeFermion(tmp_d);
  linop_single.threadedFreeFermion(tmp_s);
  for (int shift = 0; shift < nshift; shift++)
    linop_single.threadedFreeFermion(sol_s[shift]);
  
  linop.threadedFreeFermion(mtmp);
  linop.threadedFreeFermion(mmtmp);
  linop.thread_barrier();
  ////////////////////////
  return converged;
  
}

