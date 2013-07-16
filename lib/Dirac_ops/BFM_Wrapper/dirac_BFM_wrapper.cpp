/*!
 * @file dirac_BFM_wrapper.cpp
 * @brief Defines the wrapper classs for P. Boyle Bagel/BFM libs
 * Time-stamp: <2013-07-16 10:26:23 cossu>
 */

#include "dirac_BFM_wrapper.hpp"
#include "include/timings.hpp"
#include "include/messages_macros.hpp"
#include <mcheck.h>


Dirac_BFM_Wrapper_params::Dirac_BFM_Wrapper_params(XML::node BFMnode){
  XML::node top_BFM_node = BFMnode;

  XML::descend(top_BFM_node, "KernelName", MANDATORY);
  // Now specific for the scaled shamir kernel
  XML::read(top_BFM_node, "mass", mq_, MANDATORY);
  XML::read(top_BFM_node, "M5", M5_, MANDATORY);
  XML::read(top_BFM_node, "N5d", Ls_, MANDATORY);
  XML::read(top_BFM_node, "scale", scale_, MANDATORY);
}

Dirac_BFM_Wrapper_params::set_SolverParams(XML::node BFMnode){
  XML::descend(BFMnode, "Solver", MANDATORY);
  XML::read(BFMnode, "MaxIter", max_iter_, MANDATORY);
  XML::read(BFMnode, "Residual", target_, MANDATORY);
}

double Dirac_BFM_Wrapper::getMass() const{
  return BFMparams.mq_;
}

const Field* Dirac_BFM_Wrapper::getGaugeField_ptr()const{
  return u_;
}

Dirac_BFM_Wrapper::Dirac_BFM_Wrapper(XML::node node, const Field* u, DiracWilsonLike* Dirac )
  :u_(u),
   BFMparams(node),
   BFM_interface(linop),
   has_operator(false),
   has_solver_params(false),
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
  parameters.verbose = -1;
  parameters.time_report_iter=-100;

  // Threading control parameters
  threads = 64;
  bfmarg::Threads(threads);

  // Set up the diagonal part
  bfmarg::UseCGdiagonalMee(1);
  
  // Choose checkerboard
  parameters.rb_precondition_cb=Even;

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
  for(int cb=0;cb<2;cb++){
    linop.freeFermion(psi_h[cb]);
    linop.freeFermion(chi_h[cb]);
  } 
  linop.freeFermion(tmp);
}

const Field Dirac_BFM_Wrapper::mult(const Field& Fin)const{
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
  
}

Field Dirac_BFM_Wrapper::mult_test(const Field& Fin){
  if(is_initialized){
    CCIO::cout << "Dirac_BFM_Wrapper::mult_test Gauge Export\n";
    BFM_interface.GaugeExport_to_BFM(u_);
    int donrm = 0;
    int cb = Even;
    FermionField FermF(Fin);
    
    LoadSource(FermF,cb);
    CCIO::cout << "Dirac_BFM_Wrapper::mult Exporting Fermions\n";
    linop.comm_init();
#pragma omp parallel for
    for (int t=0;t<threads;t++){
      linop.MprecTilde(psi_h[cb],chi_h[cb],tmp,0,donrm);// dag=0 cb=0;
    }
    linop.comm_end();
    CCIO::cout << "Dirac_BFM_Wrapper::mult Importing Fermions\n";
    BFM_interface.FermionImport_from_BFM(FermF, chi_h[cb], cb);
    
    return FermF.data;
   
  } else {
    CCIO::cout << "The operator was not initialized yet\n";
  }
  
}



const Field Dirac_BFM_Wrapper::mult_dag(const Field& Fin)const{
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
  
}

const Field Dirac_BFM_Wrapper::md_force(const Field& eta,const Field& zeta)const{
  if(is_initialized){
    return Internal_EO->md_force(eta,zeta);
  } else {
    CCIO::cout << "The operator was not initialized yet\n";
  }

}
  /*
    Dirac_BFM_Wrapper * pThis = const_cast< Dirac_BFM_Wrapper *>(this);
    Matrix_t force[2];
    FermionField etaF(eta);
    FermionField zetaF(zeta);

    force[0] = pThis->linop.allocMatrix();
    force[1] = pThis->linop.allocMatrix();

    CCIO::cout << "Dirac_BFM_Wrapper::md_force Allocating Fermions\n";
    //Fermion_t X  = pThis->linop.threadedAllocFermion();
    //Fermion_t Y  = pThis->linop.threadedAllocFermion();
    Fermion_t X  = pThis->linop.allocFermion();
    Fermion_t Y  = pThis->linop.allocFermion();
    int cb = Even;

    CCIO::cout << "Dirac_BFM_Wrapper::md_force Gauge Export\n";
    pThis->BFM_interface.GaugeExport_to_BFM(u_);

    CCIO::cout << "Norms etaF zetaF "<< etaF.norm() << "  "<< zetaF.norm() << "\n";
    // this level of export is not needed here - change it
    CCIO::cout << "Dirac_BFM_Wrapper::md_force Exporting Fermions\n";
 
    for (int s = 0; s< BFMparams.Ls_; s++){
      pThis->BFM_interface.FermionExport_to_BFM_5D(etaF,X,cb,s); 
      pThis->BFM_interface.FermionExport_to_BFM_5D(zetaF,Y,cb,s); 
    }

    int dagno=0;
    int dagyes=1;
    CCIO::cout << "Dirac_BFM_Wrapper::md_force Barrier\n";
    //int me = pThis->linop.thread_barrier();
    pThis->linop.comm_init();
   
    CCIO::cout << "Dirac_BFM_Wrapper::md_force ZeroMatrix\n";

    
#pragma omp parallel for 
    for (int t=0;t<threads;t++){
      pThis->linop.zeroMatrix(force[0]);
      pThis->linop.zeroMatrix(force[1]);
      
      CCIO::cout << "Norms: "<< pThis->linop.norm(X) << 
	"  " << pThis->linop.norm(Y) <<"\n";
      pThis->linop.MprecDeriv(Y,X,force,dagno ,-1.0);
      //CCIO::cout << "Dirac_BFM_Wrapper::md_force MprecDeriv-2\n";
      //pThis->linop.MprecDeriv(X,Y,force,dagyes,-1.0); // no because we are taking antihermite part later
    }
    pThis->linop.comm_end();
    
    // need a gauge import
    GaugeField G_Force;
    for (int dir = 0; dir < CommonPrms::instance()->Nd(); dir++){
      pThis->BFM_interface.GaugeImport_from_BFM(&G_Force.data, force[Even], dir, Even);
      pThis->BFM_interface.GaugeImport_from_BFM(&G_Force.data, force[Odd] , dir, Odd);
    }
    G_Force.data *= 0.5;
    return G_Force.data;
    */
 
  



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
    LoadSource(source, cb);
    LoadGuess(solution, cb);

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
    GetSolution(solution,cb);
  } else {
    CCIO::cout << "The operator was not initialized yet\n";
  }
}

void Dirac_BFM_Wrapper::solve_CGNE_multishift(std::vector < FermionField > & solutions, 
					      FermionField& source,
					      vector_double shifts,
					      vector_double mresiduals){
  if(is_initialized){
 
    long double export_timing;
    FINE_TIMING_START(export_timing); 
    BFM_interface.GaugeExport_to_BFM(u_);
    FINE_TIMING_END(export_timing);
    _Message(TIMING_VERB_LEVEL, "[Timing] - Dirac_BFM_Wrapper::solve_CGNE"
	     << " - Gauge Export to BFM timing = "
	     << export_timing << std::endl); 
    int cb = Even;// by default

    std::vector< Fermion_t > ms_chi(shifts.size());
      
    vector_double m_alpha(shifts.size());
    int do_sum = 0;
    // Allocate the solutions fields in BFM
    for (int shift = 0; shift < shifts.size(); shift++){
      solutions[shift].data.resize(source.size());
      ms_chi[shift] = linop.allocFermion();//half vector
      m_alpha[shift] = 1.0;
    }
     // Load the fermion field to BFM
    LoadSource(source, cb);

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
	linop.CGNE_prec_MdagM_multi_shift(&ms_chi[0],
					  psi_h[cb],
					  &shifts[0], 
					  &m_alpha[0],
					  shifts.size(), 
					  &mresiduals[0],
					  do_sum );
      }
    }
    ///////////////////////////////////////////////
    // close communications to free the buffers
    linop.comm_end(); 
    
    // Get the solution from BFM 
    GetMultishiftSolutions(solutions,ms_chi,cb);
    
    for (int shift = 0; shift < shifts.size(); shift++){
      linop.freeFermion(ms_chi[shift]);
    }   
    
  } else {
    CCIO::cout << "The operator was not initialized yet\n";
  }
}

void Dirac_BFM_Wrapper::LoadSource(FermionField& In,int cb){
 int Ls = In.Nvol()/Nvol_;
  for (int s = 0; s< Ls; s++)
    BFM_interface.FermionExport_to_BFM_5D(In,psi_h[cb],cb,s);//valid for 4d & 5d
}

void Dirac_BFM_Wrapper::LoadGuess(FermionField& In,int cb){
 int Ls = In.Nvol()/Nvol_;
  for (int s = 0; s< Ls; s++)
    BFM_interface.FermionExport_to_BFM_5D(In,chi_h[cb],cb,s);//valid for 4d & 5d
}

void Dirac_BFM_Wrapper::GetSolution(FermionField& Out,int cb){
 int Ls = Out.Nvol()/Nvol_;
  for (int s = 0; s< Ls; s++)
    BFM_interface.FermionImport_from_BFM_5D(Out,chi_h[cb],cb,s);//valid for 4d & 5d
}

void Dirac_BFM_Wrapper::GetMultishiftSolutions(std::vector < FermionField >& Out,
					       std::vector < Fermion_t >& BFMsol,
					       int cb){
  int Ls = Out[0].Nvol()/Nvol_;
  for (int shift= 0; shift < Out.size(); shift++) 
    for (int s = 0; s< Ls; s++)
      BFM_interface.FermionImport_from_BFM_5D(Out[shift],BFMsol[shift],cb,s);//valid for 4d & 5d
  
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
  if(has_operator && has_solver_params){
    linop.init(parameters);
    //linop.comm_end();
    AllocateFields();
    is_initialized = true;
  }
}

void Dirac_BFM_Wrapper::set_SolverParams(XML::node node){
  BFMparams.set_SolverParams(node);
  parameters.max_iter = BFMparams.max_iter_;
  parameters.residual = BFMparams.target_;
  has_solver_params = true;
}
