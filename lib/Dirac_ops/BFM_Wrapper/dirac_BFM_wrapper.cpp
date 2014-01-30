/*!
 * @file dirac_BFM_wrapper.cpp
 * @brief Defines the wrapper classs for P. Boyle Bagel/BFM libs
 * Time-stamp: <2014-01-30 14:43:32 neo>
 */

#include "dirac_BFM_wrapper.hpp"

#include "include/timings.hpp"
#include "include/messages_macros.hpp"
#include <stdlib.h>     /* atoi */
#include <stdio.h>


Dirac_BFM_Wrapper_params::Dirac_BFM_Wrapper_params(XML::node BFMnode){
  XML::node top_BFM_node = BFMnode;

  XML::descend(top_BFM_node, "BFMKernel", MANDATORY);
  // Now specific for the scaled shamir kernel
  XML::read(top_BFM_node, "mass", mq_, MANDATORY);
  XML::read(top_BFM_node, "M5", M5_, MANDATORY);
  XML::read(top_BFM_node, "N5d", Ls_, MANDATORY);
  XML::read(top_BFM_node, "scale", scale_, MANDATORY);

}

Dirac_BFM_Wrapper_params::set_SolverParams(XML::node BFMnode){
  XML::read(BFMnode, "MaxIter", max_iter_, MANDATORY);
  XML::read(BFMnode, "Precision", target_, MANDATORY);
}

double Dirac_BFM_Wrapper::getMass() const{
  return BFMparams.mq_;
}

int Dirac_BFM_Wrapper::getN5() const{
  return BFMparams.Ls_;
}
const Field* Dirac_BFM_Wrapper::getGaugeField_ptr()const{
  return u_;
}

Dirac_BFM_Wrapper::Dirac_BFM_Wrapper(XML::node node, 
				     const Field* u, 
				     DiracWilsonLike_EvenOdd* Dirac,
				     Type_5d_DWF Type)
  :u_(u),
   BFMparams(node),
   BFM_interface(linop),
   has_operator(false),
   has_solver_params(false),
   is_initialized(false),
   Internal_EO(Dirac){

  // Lattice local volume
  parameters.node_latt[0]  = CommonPrms::instance()->Nx();
  parameters.node_latt[1]  = CommonPrms::instance()->Ny();
  parameters.node_latt[2]  = CommonPrms::instance()->Nz();
  parameters.node_latt[3]  = CommonPrms::instance()->Nt();

  for(int mu=0;mu<4;mu++){
    parameters.neighbour_plus[mu]  = Communicator::instance()->node_up(mu);
    parameters.neighbour_minus[mu] = Communicator::instance()->node_dn(mu);
  }

  Nvol_ =  CommonPrms::instance()->Nvol();

  // Now check if we are using multiple nodes
  // for the communications
  for(int mu=0;mu<4;mu++){
    if ( (CommonPrms::instance()->NPE(mu))>1 ) {
      parameters.local_comm[mu] = 0;
    } else {
      parameters.local_comm[mu] = 1;
    }
  }

  // Set Verbosity controls
  parameters.verbose = 0;
  parameters.time_report_iter=-100;

  // Threading control parameters
  threads = atoi(getenv("OMP_NUM_THREADS"));
  bfmarg::Threads(threads);

  // Set up the diagonal part
  bfmarg::UseCGdiagonalMee(1);
  
  // Choose checkerboard
  parameters.rb_precondition_cb=Even;

  // Check if PauliVillars
  if (Type == PauliVillars5D) BFMparams.mq_ = 1.0;

  // Set the operator depending on parameters
  set_ScaledShamirCayleyTanh(BFMparams.mq_, BFMparams.M5_, 
			      BFMparams.Ls_, BFMparams.scale_);

  // Use here a switch for assigning the
  // pointer to the solver function
  // depending on the string SolverName

  // can be changed on the fly
  //...................

  

}

Dirac_BFM_Wrapper::~Dirac_BFM_Wrapper(){
  CCIO::cout << "Destroying Dirac_BFM_Wrapper\n";
  if (is_initialized) {
    for(int cb=0;cb<2;cb++){
      linop.freeFermion(psi_h[cb]);
      linop.freeFermion(chi_h[cb]);
    } 
    linop.freeFermion(tmp);
  }
 
}

const Field Dirac_BFM_Wrapper::mult(const Field& Fin)const{
  // Now, july 2013, just a safe implementation
  return Internal_EO->mult(Fin);

  // Here comes the BFM stuff
  /*
  if(is_initialized){
    Dirac_BFM_Wrapper * pThis = const_cast< Dirac_BFM_Wrapper *>(this);
    CCIO::cout << "Dirac_BFM_Wrapper::mult Gauge Export\n";
    pThis->BFM_interface.GaugeExport_to_BFM(u_);
    int donrm = 0;
    int cb = Even;
    FermionField FermF(Fin);
    CCIO::cout << "Dirac_BFM_Wrapper::mult Exporting Fermions\n";
    pThis->BFM_interface.FermionExport_to_BFM(FermF,psi_h[cb],cb);
    pThis->linop.comm_init();
#pragma omp parallel for
    for (int t=0;t<threads;t++){
      pThis->linop.MprecTilde(psi_h[cb],chi_h[cb],tmp,0,donrm);// dag=0 cb=0;
    }
    pThis->linop.comm_end();
    CCIO::cout << "Dirac_BFM_Wrapper::mult Importing Fermions\n";
    pThis->BFM_interface.FermionImport_from_BFM(FermF, chi_h[cb], cb);
    return FermF.data;
  } else {
    CCIO::cout << "The operator was not initialized yet\n";
  }
  */
}

Field Dirac_BFM_Wrapper::mult_unprec(const Field& Fin){
  if(is_initialized){
    BFM_interface.GaugeExport_to_BFM(u_);
    int donrm = 0;
    int cb = Even;
    FermionField FermF(Fin);
    LoadSource(FermF,Even);
    LoadSource(FermF,Odd);
    linop.comm_init();
#pragma omp parallel for
    for (int t=0;t<threads;t++){
      linop.Munprec(psi_h,chi_h,tmp,0);// dag=0 
    }
    linop.comm_end();
    GetSolution(FermF,Even);
    GetSolution(FermF,Odd);
    return FermF.data;
  } else {
    CCIO::cout << "The operator was not initialized yet\n";
  }
  
}

Field Dirac_BFM_Wrapper::mult_inv_4d(const Field& Fin){
  double fact = -1.0;

  if(is_initialized){
    BFM_interface.GaugeExport_to_BFM(u_);
    int donrm = 0;
    int cb = Even;
    int dag = 0;
    FermionField FermF(Fin);
    LoadSource(FermF,Even);
    LoadSource(FermF,Odd);
    linop.comm_init();
#pragma omp parallel for
    for (int t=0;t<threads;t++){
      linop.MooeeInv(psi_h[1-cb],chi_h[1-cb],dag);// chi[odd] = Moo^-1 psi[odd]   (odd->odd) needs later
      linop.Meo(chi_h[1-cb],tmp,cb,dag);// tmp = Meo chi[odd]    (odd->even) 
      linop.axpy_norm(tmp, tmp, psi_h[cb], fact);// tmp = -1*tmp + psi[even] 
      // above -> chi[even] = psi[even] - Meo Moo^-1 psi[odd] 
      linop.MooeeInv(tmp,chi_h[cb],dag);// chi[even] = Mee^-1 tmp 
      // = Mee^-1 ( psi[even] - Meo Moo^-1 psi[odd] ) = Mee^-1 psi[even] - Mee^-1 Meo Moo^-1 psi[odd]
      
      linop.MprecTilde(chi_h[cb],psi_h[cb],tmp,1,donrm);// dag=1 
      // psi[even] = Mtilde^dag chi[even]     to get the correct inversion of Mtilde 
      linop.CGNE_prec(chi_h[cb],psi_h[cb]);// chi[even] = Mtilde^-1 psi[even]
      
      linop.Meo(chi_h[cb],tmp,1-cb,dag);// tmp = Meo chi[even]    (even->odd) 
      linop.axpy_norm(psi_h[1-cb], tmp, psi_h[1-cb], fact);// psi[odd] = -1*tmp + psi[odd]
      // = psi[odd] - Meo chi[even]
      linop.MooeeInv(psi_h[1-cb], chi_h[1-cb],dag);//chi[odd] = Moo^-1 psi[odd]
      
    }
    linop.comm_end();
    GetSolution(FermF,Even);
    GetSolution(FermF,Odd);
    return FermF.data;
  } else {
    CCIO::cout << "The operator was not initialized yet\n";
  }
  
}


Fermion_t* Dirac_BFM_Wrapper::mult_unprec_base(Fermion_t psi_in[2]){
  if(is_initialized){
    int donrm = 0;
    int cb = Even;
    BFM_interface.GaugeExport_to_BFM(u_);
    linop.comm_init();
#pragma omp parallel for
    for (int t=0;t<threads;t++){
      // Copy to the internal allocated fermion
      // to avoid cross-talking between the pointers 
      // of the PauliVillars and the DWF kernel
      // when used for the 4d BWF operator
      linop.copy(psi_h[0], psi_in[0]);
      linop.copy(psi_h[1], psi_in[1]);
      linop.Munprec(psi_h,chi_h,tmp,0);// dag=0 
    }
    linop.comm_end();
    return chi_h;
  } else {
    CCIO::cout << "The operator was not initialized yet\n";
  }
  
}

void Dirac_BFM_Wrapper::GaugeExportBFM(){
  BFM_interface.GaugeExport_to_BFM(u_);
}


Fermion_t* Dirac_BFM_Wrapper::mult_inv_4d_base(Fermion_t psi_in[2]){
  double fact = -1.0;

  if(is_initialized){
    BFM_interface.GaugeExport_to_BFM(u_);
    int donrm = 0;
    int cb = Even;
    int dag = 0;
    linop.comm_init();
#pragma omp parallel for
    for (int t=0;t<threads;t++){
     // Copy to the internal allocated fermion
      // to avoid cross-talking between the pointers 
      // of the PauliVillars and the DWF kernel
      // when used for the 4d BWF operator
      linop.copy(psi_h[0], psi_in[0]);
      linop.copy(psi_h[1], psi_in[1]);

      linop.MooeeInv(psi_h[1-cb],chi_h[1-cb],dag);// chi[odd] = Moo^-1 psi[odd]   (odd->odd) needs later
      linop.Meo(chi_h[1-cb],tmp,cb,dag);// tmp = Meo chi[odd]    (odd->even) 
      linop.axpy_norm(tmp, tmp, psi_h[cb], fact);// tmp = -1*tmp + psi[even] 
      // above -> chi[even] = psi[even] - Meo Moo^-1 psi[odd] 
      linop.MooeeInv(tmp,chi_h[cb],dag);// chi[even] = Mee^-1 tmp 
      // = Mee^-1 ( psi[even] - Meo Moo^-1 psi[odd] ) = Mee^-1 psi[even] - Mee^-1 Meo Moo^-1 psi[odd]
      
      linop.MprecTilde(chi_h[cb],psi_h[cb],tmp,1,donrm);// dag=1 
      // psi[even] = Mtilde^dag chi[even]     to get the correct inversion of Mtilde 
      linop.CGNE_prec(chi_h[cb],psi_h[cb]);// chi[even] = Mtilde^-1 psi[even]
      
      linop.Meo(chi_h[cb],tmp,1-cb,dag);// tmp = Meo chi[even]    (even->odd) 
      linop.axpy_norm(psi_h[1-cb], tmp, psi_h[1-cb], fact);// psi[odd] = -1*tmp + psi[odd]
      // = psi[odd] - Meo chi[even]
      linop.MooeeInv(psi_h[1-cb], chi_h[1-cb],dag);//chi[odd] = Moo^-1 psi[odd]
      
    }
    linop.comm_end();
    return chi_h;
  } else {
    CCIO::cout << "The operator was not initialized yet\n";
  }
  
}





const Field Dirac_BFM_Wrapper::mult_dag(const Field& Fin)const{
 // Now, july 2013, just a safe implementation
  return Internal_EO->mult_dag(Fin);


  // Here comes the BFM stuff
  /*
  if(is_initialized){
    Dirac_BFM_Wrapper * pThis = const_cast< Dirac_BFM_Wrapper *>(this);
    pThis->BFM_interface.GaugeExport_to_BFM(u_);
    int donrm = 0;
    int cb = Even;
    FermionField FermF(Fin);
    pThis->BFM_interface.FermionExport_to_BFM(FermF,psi_h[cb],cb);
    pThis->linop.comm_init();
#pragma omp parallel for
    for (int t=0;t<threads;t++){
      pThis->linop.MprecTilde(psi_h[cb],chi_h[cb],tmp,1,donrm);// dag=1 cb=0;
    }
    pThis->linop.comm_end();
    pThis->BFM_interface.FermionImport_from_BFM(FermF, chi_h[cb], cb);
    return FermF.data;
  } else {
    CCIO::cout << "The operator was not initialized yet\n";
  }
  */
}

const Field Dirac_BFM_Wrapper::md_force(const Field& eta,const Field& zeta)const{
  if(is_initialized){
    return Internal_EO->md_force(eta,zeta);
  } else {
    CCIO::cout << "The operator was not initialized yet\n";
  }

}

void Dirac_BFM_Wrapper::solve_CGNE(FermionField& solution, FermionField& source){
  if(is_initialized){
    long double export_timing;
    FINE_TIMING_START(export_timing); 
    BFM_interface.GaugeExport_to_BFM(u_);
    FINE_TIMING_END(export_timing);
    _Message(TIMING_VERB_LEVEL, "[Timing] - Dirac_BFM_Wrapper::solve_CGNE"
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

    // need to initialize the comms because of the buffer 
    // assignments by BGNET library
    // the order of comm_init can broke previous initializations
    // so we are doing right now.
    linop.comm_init();

    //////////////////////////// Execute the solver
    #pragma omp parallel
    {
    #pragma omp for 
      for (int t=0;t<threads;t++){
	linop.CGNE_prec(chi_h[cb],psi_h[cb]);
      }
    }
    ///////////////////////////////////////////////

    // close communications to free the buffers
    linop.comm_end();
    
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
void Dirac_BFM_Wrapper::solve_CGNE_multishift(std::vector < FermionField > & solutions, 
                                              FermionField& source,
                                              vector_double shifts,
                                              vector_double mresiduals){
  if(is_initialized){
    std::vector < FermionField > internal_sol(shifts.size());
    FermionField BFMsource(Nvol_*BFMparams.Ls_);
    Fermion_t ms_chi[shifts.size()];
    vector_double m_alpha(shifts.size());
    long double export_timing;
    int cb = Even;// by default
    int dont_sum = 0;
    int half_vec = source.size()/BFMparams.Ls_;

    // Import gauge field
    BFM_interface.GaugeExport_to_BFM(u_);

    // Allocate the solutions fields in BFM
    for (int shift = 0; shift < shifts.size(); shift++){
      solutions[shift].data.resize(source.size());     //half vector
      internal_sol[shift].data.resize(2*source.size());// full vector
      ms_chi[shift] = linop.allocFermion();            //half vector
      m_alpha[shift] = 1.0; //fixed
    }

    // assumes even/odd indexing
    for (int s =0 ; s< BFMparams.Ls_; s++){
      
	memcpy((void*)&(BFMsource.data[2*s*half_vec]),
	       (void*)&(source.data[s*half_vec]), 
	       sizeof(double)*half_vec);
      
	/*	
      for (int i = 0; i < half_vec; i++){
	BFMsource.data.set(i+    2*s*half_vec, source.data[i+half_vec*s]);//just even part is enough
	BFMsource.data.set(i+(2*s+1)*half_vec, 0.0);//just even part is enough
	}*/
	
    }
    // Load the fermion field to BFM
    LoadSource(BFMsource, cb);

    // need to initialize the comms because of the buffer 
    // assignments by BGNET library
    // the order of comm_init can broke previous initializations
    // so we are doing right now.
    linop.comm_init();

   //////////////////////////// Execute the solver
    #pragma omp parallel for
    for (int t=0;t<threads;t++){
      int iter = linop.CGNE_prec_MdagM_multi_shift(ms_chi,
					psi_h[cb],
					&shifts[0], 
					&m_alpha[0],
					shifts.size(), 
					&mresiduals[0],
					dont_sum);
      
    }
    
    ///////////////////////////////////////////////
    // close communications to free the buffers
    linop.comm_end(); 
     
    // Get the solution from BFM and convert basis 
    GetMultishiftSolutions(internal_sol,ms_chi,cb);
      
    for (int shift = 0; shift < shifts.size(); shift++){
      for (int s =0 ; s< BFMparams.Ls_; s++){
	
	memcpy((void*)&(solutions[shift].data[half_vec*s]),
	       (void*)&(internal_sol[shift].data[2*s*half_vec]), 
	       sizeof(double)*half_vec);
	/*
	for (int i = 0; i < half_vec; i++){
	  solutions[shift].data.set(i+half_vec*s, internal_sol[shift].data[i+2*s*half_vec]);//copy only the even part
	}
	*/
      }
      linop.freeFermion(ms_chi[shift]);
    }
    
  } else {
    CCIO::cout << "The operator was not initialized yet\n";
  }
}

void Dirac_BFM_Wrapper::LoadSource(FermionField& In,int cb){
  for (int s = 0; s< BFMparams.Ls_; s++)
    BFM_interface.FermionExport_to_BFM_5D(In,psi_h[cb],cb,s);//valid for 4d & 5d
}

Fermion_t* Dirac_BFM_Wrapper::LoadFullSource(FermionField& In){
  for (int s = 0; s< BFMparams.Ls_; s++){
    BFM_interface.FermionExport_to_BFM_5D(In,psi_h[Even],Even,s);//valid for 4d & 5d
    BFM_interface.FermionExport_to_BFM_5D(In,psi_h[Odd],Odd,s);//valid for 4d & 5d
  }
  return psi_h;
}

void Dirac_BFM_Wrapper::LoadGuess(FermionField& In,int cb){
  for (int s = 0; s< BFMparams.Ls_; s++)
    BFM_interface.FermionExport_to_BFM_5D(In,chi_h[cb],cb,s);//valid for 4d & 5d
}

void Dirac_BFM_Wrapper::GetSolution(FermionField& Out,int cb){
 int Ls = BFMparams.Ls_;
  for (int s = 0; s< Ls; s++)
    BFM_interface.FermionImport_from_BFM_5D(Out,chi_h[cb],cb,s, Ls);//valid for 4d & 5d
}

void Dirac_BFM_Wrapper::GetFullSolution(FermionField& Out){
 int Ls = BFMparams.Ls_;
 for (int s = 0; s< Ls; s++){
    BFM_interface.FermionImport_from_BFM_5D(Out,chi_h[Even],Even,s, Ls);//valid for 4d & 5d
    BFM_interface.FermionImport_from_BFM_5D(Out,chi_h[Odd] ,Odd ,s, Ls);//valid for 4d & 5d
 }
}


void Dirac_BFM_Wrapper::GetMultishiftSolutions(std::vector < FermionField >& Out,
					       Fermion_t *BFMsol,
					       int cb){
   for (int shift= 0; shift < Out.size(); shift++) 
    for (int s = 0; s< BFMparams.Ls_; s++)
      BFM_interface.FermionImport_from_BFM_5D(Out[shift],
					      BFMsol[shift],
					      cb,s, BFMparams.Ls_);//valid for 4d & 5d
   
}



void  Dirac_BFM_Wrapper::AllocateFields(){
  if(has_operator){
    // Allocate the fermions
    for(int cb=0;cb<2;cb++){
      psi_h[cb] = linop.allocFermion();//half vector
      CCIO::cout << "BFM phi_h["<<cb<<"] allocated \n";
      chi_h[cb] = linop.allocFermion();
      CCIO::cout << "BFM chi_h["<<cb<<"] allocated \n";
    }
    tmp = linop.allocFermion();// temporary
  }
}



//////////////////////////////// BFM Operators
void Dirac_BFM_Wrapper::set_ScaledShamirCayleyTanh(double mq, 
						   double M5, 
						   int Ls, 
						   double scale){
  // Note M5 has different sign from our convention
  parameters.ScaledShamirCayleyTanh(mq,-M5,Ls,scale);
  has_operator = true;
}

void Dirac_BFM_Wrapper::initialize(){
  // Initialize the linear operator 
  if(has_operator && has_solver_params && !is_initialized){
    linop.init(parameters);
    AllocateFields();
    is_initialized = true;
  }
}

void Dirac_BFM_Wrapper::set_SolverParams(XML::node node){
  BFMparams.set_SolverParams(node);
  parameters.max_iter = BFMparams.max_iter_;
  parameters.residual = sqrt(BFMparams.target_);
  has_solver_params = true;
}


size_t Dirac_BFM_Wrapper::fsize() const {
  return Internal_EO->fsize();
}

size_t Dirac_BFM_Wrapper::gsize() const {
  return Internal_EO->gsize();
}
