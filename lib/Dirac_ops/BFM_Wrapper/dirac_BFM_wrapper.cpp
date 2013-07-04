/*!
 * @file dirac_BFM_wrapper.cpp
 * @brief Defines the wrapper classs for P. Boyle Bagel/BFM libs
 * Time-stamp: <2013-07-03 18:02:31 cossu>
 */

#include "dirac_BFM_wrapper.hpp"


Dirac_BFM_Wrapper::Dirac_BFM_Wrapper(const Field* u)
  :u_(u),BFM_interface(linop),
   has_operator(false),
   has_solver_params(false){
  //Initializes the BFM classes
  
  // Lattice local volume
  parameters.node_latt[0]  = CommonPrms::instance()->Nx();
  parameters.node_latt[1]  = CommonPrms::instance()->Ny();
  parameters.node_latt[2]  = CommonPrms::instance()->Nz();
  parameters.node_latt[3]  = CommonPrms::instance()->Nt();
  Nvol_ =  CommonPrms::instance()->Nvol();
  // Checks if we are using multiple nodes
  // for the communications
  for(int mu=0;mu<4;mu++){
    if ( (CommonPrms::instance()->NPE(mu))>1 ) {
      parameters.local_comm[mu] = 0;
    } else {
      parameters.local_comm[mu] = 1;
    }
  }

  // Verbosity controls
  parameters.verbose = 0;
  parameters.time_report_iter=-100;

  // Threading control
  threads = 64;
  bfmarg::Threads(threads);

  // Sets up the diagonal part
  bfmarg::UseCGdiagonalMee(1);

  // Checkerboard
  parameters.rb_precondition_cb=Even;

  // Use here a switch for assigning the
  // pointer to the solver function
  // depending on the string SolverName

  // can be changed on the fly

}

Dirac_BFM_Wrapper::~Dirac_BFM_Wrapper(){
  for(int cb=0;cb<2;cb++){
    linop.freeFermion(psi_h[cb]);
    linop.freeFermion(chi_h[cb]);
  } 
  linop.freeFermion(tmp);
}

const Field Dirac_BFM_Wrapper::mult(const Field& Fin)const{
  if(has_operator){
    Dirac_BFM_Wrapper * pThis = const_cast< Dirac_BFM_Wrapper *>(this);
      pThis->BFM_interface.GaugeExport_to_BFM(u_);
    int donrm = 0;
    int cb = Even;
    FermionField FermF(Fin);
    pThis->BFM_interface.FermionExport_to_BFM(FermF,psi_h[cb],cb);
    
#pragma omp parallel for
    for (int t=0;t<threads;t++){
      pThis->linop.MprecTilde(psi_h[cb],chi_h[cb],tmp,0,donrm);// dag=0 cb=0;
    }
    
    pThis->BFM_interface.FermionImport_from_BFM(FermF, chi_h[cb], cb);
    return FermF.data;
  } else {
    CCIO::cout << "The operator was not set\n";
  }
  
}

const Field Dirac_BFM_Wrapper::mult_dag(const Field& Fin)const{
  if(has_operator){
    Dirac_BFM_Wrapper * pThis = const_cast< Dirac_BFM_Wrapper *>(this);
    pThis->BFM_interface.GaugeExport_to_BFM(u_);
    int donrm = 0;
    int cb = Even;
    FermionField FermF(Fin);
    pThis->BFM_interface.FermionExport_to_BFM(FermF,psi_h[cb],cb);
#pragma omp parallel for
    for (int t=0;t<threads;t++){
      pThis->linop.MprecTilde(psi_h[cb],chi_h[cb],tmp,1,donrm);// dag=1 cb=0;
    }
    
    pThis->BFM_interface.FermionImport_from_BFM(FermF, chi_h[cb], cb);
    return FermF.data;
  } else {
    CCIO::cout << "The operator or the solver parameters are not set\n";
  }
  
}

void Dirac_BFM_Wrapper::solve_CGNE(FermionField& solution, FermionField& source){
  if(has_operator && has_solver_params){
    BFM_interface.GaugeExport_to_BFM(u_);
    int cb = Even;// by default
    LoadField(source, cb);
    #pragma omp parallel
    {
    #pragma omp for 
      for (int t=0;t<threads;t++){
	linop.CGNE_prec(chi_h[cb],psi_h[cb]);
      }
    }
    GetField(solution,cb);
  } else {
    CCIO::cout << "The operator is not set\n";
  }
}

void Dirac_BFM_Wrapper::LoadField(FermionField& In,int cb){
 int Ls = In.Nvol()/Nvol_;
  for (int s = 0; s< Ls; s++)
    BFM_interface.FermionExport_to_BFM_5D(In,psi_h[cb],cb,s);//valid for 4d & 5d
}

void Dirac_BFM_Wrapper::GetField(FermionField& Out,int cb){
 int Ls = Out.Nvol()/Nvol_;
  for (int s = 0; s< Ls; s++)
    BFM_interface.FermionImport_from_BFM_5D(Out,chi_h[cb],cb,s);//valid for 4d & 5d
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
void Dirac_BFM_Wrapper::set_ScaledShamirCayleyTanh(double mq, double M5, int Ls, double scale){
  // Note M5 has different sign from our convention
  parameters.ScaledShamirCayleyTanh(mq,-M5,Ls,scale);
  has_operator = true;
}

void Dirac_BFM_Wrapper::initialize(){
  // Initialize the linear operator 
  if(has_operator && has_solver_params){
    linop.init(parameters);
    AllocateFields();
  }
}

void Dirac_BFM_Wrapper::solver_params(int max, double target_residual){
  parameters.max_iter=max;
  parameters.residual=target_residual; 
  has_solver_params = true;
}
