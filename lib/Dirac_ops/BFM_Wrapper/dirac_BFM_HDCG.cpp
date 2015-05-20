/*!
 * @file dirac_BFM_HDCG_wrapper.cpp
 * @brief Defines the wrapper class methods for P. Boyle HDCG inverter
 * Time-stamp: <2015-05-20 13:57:48 cossu>
 */
#include <omp.h>
#include <math.h>
#include "dirac_BFM_HDCG.hpp"
#include "include/commonPrms.hpp"
#include "Communicator/communicator.hpp"
#include "include/timings.hpp"
#include "include/messages_macros.hpp"
#include "include/errors.hpp"

void Dirac_BFM_HDCG_Wrapper::open_comms(){
  linop_r.comm_init();
  linop_single.comm_init();
  linop.comm_init();
}


void Dirac_BFM_HDCG_Wrapper::close_comms(){
  linop.comm_end();
  linop_single.comm_end();
  linop_r.comm_end();
}


void Dirac_BFM_HDCG_Wrapper::initialize(XML::node node){
  HDCG_init(node);
  
  long double subspace_timer, compute_timer, refine_timer;
  FINE_TIMING_START(subspace_timer);
  
  // 1. Init subspace (computes the vectors)
  CCIO::cout << ".::::::::::::::  Init subspace\n";
  HDCG_subspace_init();

  FINE_TIMING_END(subspace_timer);
  _Message(TIMING_VERB_LEVEL,"Timer subspace: " << subspace_timer << " seconds \n");

  FINE_TIMING_START(compute_timer);
  // 2. Compute the subspace (coarse space Little Dirac operator construction)
  CCIO::cout << ".::::::::::::::  Compute subspace\n";
  HDCG_subspace_compute(0);  // 0 = double, 1 = single
  
  FINE_TIMING_END(compute_timer);
  _Message(TIMING_VERB_LEVEL,"Timer Little Dirac compute: " << compute_timer << " seconds \n");

  // 3. Refine the subspace (optional)
  if (do_refine()){

    CCIO::cout << ".::::::::::::::  Refine subspace\n";
    FINE_TIMING_START(refine_timer);
    HDCG_subspace_refine(); 
    HDCG_subspace_compute(0);  // 0 = double, 1 = single
    FINE_TIMING_END(refine_timer);
    _Message(TIMING_VERB_LEVEL,"Timer Little Dirac refine: " << refine_timer << " seconds \n");
  }

 
}

void Dirac_BFM_HDCG_Wrapper::HDCG_init(XML::node HDCGnode) {

  //Get the HDCG specific parameters from the XML node
  XML::node top_HDCG_node = HDCGnode;
  XML::descend(top_HDCG_node, "HDCG_Solver", MANDATORY);
  
  //////////////////////////////////////////////// Set parameters
  // Basic
  XML::read(top_HDCG_node, "NumberSubspace", HDCGParams.NumberSubspace, MANDATORY);
  HDCGParams.Ls = BFMparams.Ls_;
  std::vector<int> blockdim;
  XML::read_array(top_HDCG_node, "BlockDim", blockdim, MANDATORY);
  HDCGParams.Block[0] = blockdim[0];
  HDCGParams.Block[1] = blockdim[1];
  HDCGParams.Block[2] = blockdim[2];
  HDCGParams.Block[3] = blockdim[3];
  HDCGParams.Block[4] = blockdim[4];

  // First pass subspace controls  - Subspace Rationals
  XML::read(top_HDCG_node, "SubspaceRationalLs", HDCGParams.SubspaceRationalLs, MANDATORY);
  XML::read(top_HDCG_node, "SubspaceRationalLo", HDCGParams.SubspaceRationalLo, MANDATORY);
  XML::read(top_HDCG_node, "SubspaceRationalMass", HDCGParams.SubspaceRationalMass, MANDATORY);
  XML::read(top_HDCG_node, "SubspaceRationalResidual", HDCGParams.SubspaceRationalResidual, MANDATORY);
  HDCGParams.SubspaceSurfaceDepth = HDCGParams.SubspaceRationalLs;//default value if not specified
  XML::read(top_HDCG_node, "SubspaceSurfaceDepth", HDCGParams.SubspaceSurfaceDepth); 

  // Second pass subspace controls
  XML::read(top_HDCG_node, "SubspaceRationalRefine", HDCGParams.SubspaceRationalRefine, MANDATORY);
  HDCGParams.SubspaceRationalRefineLo =   0.001;// default value if not specified
  XML::read(top_HDCG_node, "SubspaceRationalRefineLo", HDCGParams.SubspaceRationalRefineLo);
  HDCGParams.SubspaceRationalRefineResidual =   1.0e-5;// default value if not specified
  XML::read(top_HDCG_node, "SubspaceRationalRefineResidual", HDCGParams.SubspaceRationalRefineResidual);

  // Little Dirac operator deflation controls
  XML::read(top_HDCG_node, "LdopDeflVecs", HDCGParams.LdopDeflVecs, MANDATORY);
  XML::read(top_HDCG_node, "LittleDopSolverResidualSubspace", HDCGParams.LittleDopSolverResidualSubspace, MANDATORY);
  XML::read(top_HDCG_node, "LittleDopSolverResidualInner", HDCGParams.LittleDopSolverResidualInner, MANDATORY);
  XML::read(top_HDCG_node, "LittleDopSubspaceRational", HDCGParams.LittleDopSubspaceRational, MANDATORY);
  XML::read(top_HDCG_node, "LittleDopSolverResidualVstart", HDCGParams.LittleDopSolverResidualVstart, MANDATORY);
  XML::read(top_HDCG_node, "LittleDopSolverIterMax", HDCGParams.LittleDopSolverIterMax, MANDATORY);
  //////////////////////////////////////////////////////////

  // Outer solver parameters
  XML::read(top_HDCG_node, "PreconditionerKrylovResidual", HDCGParams.PreconditionerKrylovResidual, MANDATORY);
  XML::read(top_HDCG_node, "PreconditionerKrylovIterMax", HDCGParams.PreconditionerKrylovIterMax, MANDATORY);
  XML::read(top_HDCG_node, "PreconditionerKrylovShift", HDCGParams.PreconditionerKrylovShift, MANDATORY);
  HDCGParams.PcgSingleShift = 0.0; // default value if not specified
  XML::read(top_HDCG_node, "PcgSingleShift", HDCGParams.PcgSingleShift);

  std::string LittleDopSolver_str;
  XML::read(top_HDCG_node, "LittleDopSolver", LittleDopSolver_str, MANDATORY);
  if (!LittleDopSolver_str.compare("ADef2")){
    HDCGParams.LittleDopSolver = LittleDopSolverADef2;
  }  else if (!LittleDopSolver_str.compare("DeflCG")){
    HDCGParams.LittleDopSolver = LittleDopSolverDeflCG;
  } else   if (!LittleDopSolver_str.compare("CG")){
    HDCGParams.LittleDopSolver = LittleDopSolverCG;
  } else   if (!LittleDopSolver_str.compare("MCR")){
    HDCGParams.LittleDopSolver = LittleDopSolverMCR;
  } else {
    Errors::ParameterErr("Undefined LittleDopSolver declaration. Check XML\n");
  }

  std::string LittleDopM1control_str;
  XML::read(top_HDCG_node, "LittleDopM1control", LittleDopM1control_str, MANDATORY);
  if (!LittleDopM1control_str.compare("Chebyshev")){
    HDCGParams.LdopM1control = LdopM1Chebyshev;
  }  else   if (!LittleDopM1control_str.compare("Mirs")){
    HDCGParams.LdopM1control = LdopM1Mirs;
  }  else   if (!LittleDopM1control_str.compare("MirsPoly")){
    HDCGParams.LdopM1control = LdopM1MirsPoly;
  }  else   if (!LittleDopM1control_str.compare("MirsPolyRecord")){
    HDCGParams.LdopM1control = LdopM1MirsPolyRecord;
  }  else   if (!LittleDopM1control_str.compare("DiagInv")){
    HDCGParams.LdopM1control = LdopM1DiagInv;
  }  else {
    Errors::ParameterErr("Undefined LittleDopM1control declaration. Check XML\n");
  }
  XML::read(top_HDCG_node, "LdopM1Lo", HDCGParams.LdopM1Lo, MANDATORY);
  XML::read(top_HDCG_node, "LdopM1Hi", HDCGParams.LdopM1Hi, MANDATORY);
  XML::read(top_HDCG_node, "LdopM1iter", HDCGParams.LdopM1iter, MANDATORY);
  XML::read(top_HDCG_node, "LdopM1resid", HDCGParams.LdopM1resid, MANDATORY);


  XML::read(top_HDCG_node, "Flexible", HDCGParams.Flexible, MANDATORY);
  if (HDCGParams.Flexible)
    BFMparams.is_mixed_precision = true; // forces the initialization of single precision gauge field in solvers

  // Complicated init sequence
  bfmarg bfma;
  bfmActionParams *bfmap = (bfmActionParams *) &bfma;

  // Physical params
  *bfmap = parameters; // from the parent class that has been already initialized by XML

  // Algorithm & code control
  bfma.Threads(omp_get_max_threads());
  bfma.Verbose(BfmMessage|BfmDebug|BfmPerformance|BfmError);
  bfma.time_report_iter  =-100;

  // Local geometry, get from IroIro environment
  bfma.node_latt[0] =  CommonPrms::instance()->Nx();
  bfma.node_latt[1] =  CommonPrms::instance()->Ny();
  bfma.node_latt[2] =  CommonPrms::instance()->Nz();
  bfma.node_latt[3] =  CommonPrms::instance()->Nt();

  int procs[4];
  for(int mu=0;mu<4;mu++){
    procs[mu] = CommonPrms::instance()->NPE(mu);
    if (procs[mu]>1) bfma.local_comm[mu] = 0;
    else             bfma.local_comm[mu] = 1;
  }

  // Node topology
  for(int mu=0;mu<4;mu++){
    bfma.neighbour_plus[mu]  = Communicator::instance()->node_up(mu);
    bfma.neighbour_minus[mu] = Communicator::instance()->node_dn(mu);
  }
  

  // Initialization of internal operators
  linop_single.init(bfma);
  linop_single.BossLog("HDCG wrapper inititialised single precision linop, parameters %d %f\n", bfma.Ls, bfma.mass);
 
  linop.init(bfma);
  linop.BossLog("HDCG wrapper inititialised double precision linop, parameters %d %f\n", bfma.Ls, bfma.mass);
  
  bfma.Ls   = HDCGParams.SubspaceRationalLs;
  bfma.mass = HDCGParams.SubspaceRationalMass;
  
  linop_r.BossLog("HDCG inititialising subspace generation Ls=%d mass %le\n",bfma.Ls,bfma.mass);
  linop_r.init(bfma);

  // Merge the Params into constructor.
  int GlobalSize[4];
  int SubGridSize[4];
  int NodeCoord[4];
  
  for (int mu =0; mu< 4; mu++){
    GlobalSize[mu]  = CommonPrms::instance()->global_size(mu);
    SubGridSize[mu] = CommonPrms::instance()->local_size(mu);
    NodeCoord[mu]   = Communicator::instance()->ipe(mu);
  }

  linop.BossLog("HDCG class initialization\n");
  ldop_d = new BFM_HDCG_Extend<double>(HDCGParams.Ls,HDCGParams.NumberSubspace,HDCGParams.Block,HDCGParams.Block,GlobalSize, SubGridSize, NodeCoord, &linop,&linop_single);
  ldop_d->SetParams(HDCGParams);

  //Allocate fields
  AllocateFields(); // just for the double precision fields
  is_initialized = true;
}

void Dirac_BFM_HDCG_Wrapper::HDCG_reinit(){
  is_initialized = false;
  delete(ldop_d);

  bfmarg bfma;
  bfmActionParams *bfmap = (bfmActionParams *) &bfma;

  // Physical params
  *bfmap = parameters; // from the parent class that has been already initialized by XML

  // Algorithm & code control
  bfma.Threads(omp_get_max_threads());
  bfma.Verbose(BfmMessage|BfmDebug|BfmPerformance|BfmError);
  bfma.time_report_iter  =-100;

  // Local geometry, get from IroIro environment
  bfma.node_latt[0] =  CommonPrms::instance()->Nx();
  bfma.node_latt[1] =  CommonPrms::instance()->Ny();
  bfma.node_latt[2] =  CommonPrms::instance()->Nz();
  bfma.node_latt[3] =  CommonPrms::instance()->Nt();

  int procs[4];
  for(int mu=0;mu<4;mu++){
    procs[mu] = CommonPrms::instance()->NPE(mu);
    if (procs[mu]>1) bfma.local_comm[mu] = 0;
    else             bfma.local_comm[mu] = 1;
  }

  // Node topology
  for(int mu=0;mu<4;mu++){
    bfma.neighbour_plus[mu]  = Communicator::instance()->node_up(mu);
    bfma.neighbour_minus[mu] = Communicator::instance()->node_dn(mu);
  }

  bfma.Ls   = HDCGParams.SubspaceRationalLs;
  bfma.mass = HDCGParams.SubspaceRationalMass;
  
  int GlobalSize[4];
  int SubGridSize[4];
  int NodeCoord[4];
  
  for (int mu =0; mu< 4; mu++){
    GlobalSize[mu]  = CommonPrms::instance()->global_size(mu);
    SubGridSize[mu] = CommonPrms::instance()->local_size(mu);
    NodeCoord[mu]   = Communicator::instance()->ipe(mu);
  }



  linop_r.BossLog("HDCG inititialising subspace generation Ls=%d mass %le\n",bfma.Ls,bfma.mass);
  linop_r.init(bfma);
 
  linop.BossLog("HDCG class initialization\n");
  ldop_d = new BFM_HDCG_Extend<double>(HDCGParams.Ls,HDCGParams.NumberSubspace,HDCGParams.Block,HDCGParams.Block,GlobalSize, SubGridSize, NodeCoord, &linop,&linop_single);
  ldop_d->SetParams(HDCGParams); 

  is_initialized = true;

}

void Dirac_BFM_HDCG_Wrapper::HDCG_subspace_rotate(double *phases_by_pi,double *phases_dir, int sgn){
  
  /* TwistedBC does the following
     for(int i=0; i < Nd-1; i++) {
     int direction = phases_dir[i]; 
     Real arg = phases_by_pi[i]*onepi / Real( nrow[ direction ] );
     Real re = cos(arg);
     Real im = sin(arg);
     Complex phase = cmplx( re, im );
     u[ direction ] *= phase;
     }
     * Suppose 
     *    D psi = 0
     * Want to modify subspace vector psi so that
     *    D^tw psi^tw = 0
     * If 
     *    psi^tw(x) = e^{-i x arg} psi(x)
     * Then 
     *    U_\mu^tw psi^tw(x+mu) = e^{i x arg} U_\mu psi(x+mu) 
     */
  linop.BossLog("HDCG Rotating subspace for twisted BCs\n");
  

  double onepi = M_PI;
  int LocalVol = CommonPrms::instance()->Nvol();
  int Ls=linop.Ls;  

  vector_double theta, eitheta_re, eitheta_im;// scalar fields
  theta.resize(LocalVol); 
  eitheta_re.resize(LocalVol); 
  eitheta_im.resize(LocalVol); 
  for(int i=0; i < ND_-1; i++) { 
    int mu = phases_dir[i];
    
    for(int site = 0; site < LocalVol; site++) {
	int coord[4];
	coord[0] = SiteIndex::instance()->c_x(site);
	coord[1] = SiteIndex::instance()->c_y(site);
	coord[2] = SiteIndex::instance()->c_z(site);
	coord[3] = SiteIndex::instance()->c_t(site);
	
	theta[site] = theta[site] + (SiteIndex::instance()->*SiteIndex::global_idx[mu])(coord[mu]) * phases_by_pi[i]* onepi / CommonPrms::instance()->global_size(mu) ;

	eitheta_re[site] = cos(theta[site]);
	eitheta_im[site] = -sgn*sin(theta[site]);
	
	
      }
	  }
    
    FermionField copy(LocalVol*Ls); // 5d object
  
    for(int v=0;v<ldop_d->Nvec;v++){
      BFM_interface.FermionExport_to_BFM(copy, ldop_d->subspace_d[v], 1); //last parameter Even/Odd
      //linop.exportFermion(copy,ldop_d->subspace_d[v],1);
      for(int s=0;s<Ls;s++){
	for(int site = 0; site < LocalVol; site++) {
	    for (int n_in = 0; n_in < copy.Nin()/2; n_in ++){
	      int idx = copy.format.index(n_in, site+LocalVol*s);
	      int idxp = copy.format.index(n_in+1, site+LocalVol*s);
	      copy.data.set(idx, eitheta_re[site]*copy[idx]);
	      copy.data.set(idxp, eitheta_im[site]*copy[idxp]);
	    }
	  }
	      }
	      BFM_interface.FermionImport_from_BFM(copy, ldop_d->subspace_d[v], 1); //last parameter Even/Odd     
	    //linop.importFermion(copy,ldop_d->subspace_d[v],1);
	    
	
      }
  
    
  }


void Dirac_BFM_HDCG_Wrapper::HDCG_set_mass(double mass){
    linop.BossLog("HDCG setting mass %le\n",mass);
    linop.GeneralisedFiveDimEnd();
    linop.mass = mass;
    linop.GeneralisedFiveDimInit();

    linop_single.GeneralisedFiveDimEnd();
    linop_single.mass = mass;
    linop_single.GeneralisedFiveDimInit();
}


void Dirac_BFM_HDCG_Wrapper::HDCG_subspace_init(){

  //int sloppy_relax=HDCGParams.SubspaceRationalRefine;
  int sloppy_relax=1;

  uint32_t seed = (ldop_d->ncoor[0]*17    + 
		   ldop_d->ncoor[1]*251   +
		   ldop_d->ncoor[2]*1023  + 
		   ldop_d->ncoor[3]*8191);
  
  srand48(seed);
  
  if(is_initialized){
    
    if ( sloppy_relax ) {
      CCIO::cout << "HDCG_subspace_init sloppy relax \n";
      BFM_interface.GaugeExport_to_BFM(u_);
      BFM_interface_r.GaugeExport_to_BFM(u_);
      ldop_d->RelaxSubspace<rFloat>(&linop_r);
    } else {
      CCIO::cout << "HDCG_subspace_init Relax in double \n";
      BFM_interface.GaugeExport_to_BFM(u_);
      BFM_interface_r.GaugeExport_to_BFM(u_);
      ldop_d->RelaxSubspace<double>(&linop);
    }

    CCIO::cout << "HDCG_subspace_init SinglePrecSubspace \n";
    ldop_d->SinglePrecSubspace();

  } else {
    CCIO::cout << "The operator was not initialized yet\n";
  }

}

void Dirac_BFM_HDCG_Wrapper::HDCG_subspace_refine(){
  
  int Nvec = HDCGParams.NumberSubspace;
  
  Fermion_t src = linop.allocFermion();
  Fermion_t sol = linop.allocFermion();
  Fermion_t tmp = linop.allocFermion();
  
  int threads = linop.threads;
  BFM_HDCG_Extend<float>  *ldop_f;
  int GlobalSize[4];
  int SubGridSize[4];
  int NodeCoord[4];
  
  for (int mu =0; mu< 4; mu++){
    GlobalSize[mu] = CommonPrms::instance()->global_size(mu);
    SubGridSize[mu] = CommonPrms::instance()->local_size(mu);
    NodeCoord[mu] = Communicator::instance()->ipe(mu);
  }
  ldop_f = new BFM_HDCG_Extend<float>(HDCGParams.Ls,HDCGParams.NumberSubspace,HDCGParams.Block,HDCGParams.Block,GlobalSize, SubGridSize, NodeCoord, &linop,&linop_single);
  printf("ldop_f node number %d\n",MyNodeNumber());
  ldop_f->CloneSubspace<double>(*ldop_d);
  ldop_f->SetParams(HDCGParams);

  ldop_f->PcgSingleShift = HDCGParams.SubspaceRationalRefineLo;
  linop.residual = HDCGParams.SubspaceRationalRefineResidual;
  ldop_f->PcgType        = PcgAdef2fSingleShift;
  
#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<threads;i++) {
      for(int v=0;v<Nvec;v++){
	
	double nn=linop.norm(ldop_d->subspace_r[v]);
	linop.ThreadBossLog("HDCG Refining vector[%d] of norm %le\n",v,nn);
	
	linop.copy(src,ldop_d->subspace_r[v]);
	linop.set_zero(sol);
	ldop_f->Pcg(ldop_d,sol,src,tmp);
	linop.copy(ldop_d->subspace_r[v],sol);
	
	nn=linop.norm(ldop_d->subspace_r[v]);
	//	  linop.ThreadBossLog(" Refined vector[%d] new norm %le\n",v,nn);
      }
      for(int v=0;v<Nvec;v++){
	linop.copy(ldop_d->subspace_d[v],ldop_d->subspace_r[v]);
      }
    }
  }
  ldop_f->end();
  delete ldop_f;
  
  linop.freeFermion(tmp);
  linop.freeFermion(src);
  linop.freeFermion(sol);
  
  ldop_d->FreeSingleSubspace();
  ldop_d->OrthogonaliseSubspace();
  ldop_d->SinglePrecSubspace();
  
}

void Dirac_BFM_HDCG_Wrapper::HDCG_subspace_free(void){
  ldop_d->FreeDoubleSubspace();
  ldop_d->FreeSingleSubspace();
}


void Dirac_BFM_HDCG_Wrapper::HDCG_subspace_compute(int sloppy){
  
  int threads = linop.threads;
  int Nvec = HDCGParams.NumberSubspace;
  
  linop.BossMessage("Computing subspace of %d vectors\n",HDCGParams.NumberSubspace);
  //open_comms();

  if ( sloppy ) { 
    ldop_d->ComputeLittleMatrixColored<float>(&linop_single,ldop_d->subspace_f);
    ldop_d->LdopDeflationBasisSize=0;
  } else { 
    ldop_d->ComputeLittleMatrixColored<double>(&linop,ldop_d->subspace_d);
    // NB due to BG/Q QPX have hard coded DeflationBasisSize must be multiple of 4.
    // This may become even more on future machines.
#if 1
    if( HDCGParams.SubspaceSurfaceDepth < HDCGParams.Ls/2) {
      linop.BossMessage("Zeroing s=[%d,%d]\n",HDCGParams.SubspaceSurfaceDepth,HDCGParams.Ls-HDCGParams.SubspaceSurfaceDepth-1);
#pragma omp parallel 
      {
#pragma omp for 
	for(int i=0;i<threads;i++) {
	  for(int v=0;v<Nvec;v++){
	    double nno = linop.norm(ldop_d->subspace_d[v]);
	    for(int s=HDCGParams.SubspaceSurfaceDepth; s<HDCGParams.Ls-HDCGParams.SubspaceSurfaceDepth;s++){
	      linop.axpby_ssp(ldop_d->subspace_d[v],0.0,ldop_d->subspace_d[v],0.0,ldop_d->subspace_d[v],s,s);
	    }
	    double nnr = linop.norm(ldop_d->subspace_d[v]);
	      linop.ThreadBossLog("Restricted vec[%d] %le -> %le \n",v,nno,nnr);
	  }
	}
      }
    }
#endif
    
    uint32_t seed = ldop_d->ncoor[0]*11 + ldop_d->ncoor[1]*161
      + ldop_d->ncoor[2]*357 + ldop_d->ncoor[3]*1203;
    srand48(seed);
    
    int LdopDeflVecs = HDCGParams.LdopDeflVecs;
    ldop_d->LdopDeflationBasisSize=0;
    //      ldop_d->LdopDeflationBasisTrivial(); DEBUG
    int min =4;
      if ( ldop_d->LdopDeflationBasisSize> min )  min = ldop_d->LdopDeflationBasisSize;
      ldop_d->LdopDeflationIsDiagonal = 0; // force reinitialization      
      for (int v = min ; v<=LdopDeflVecs; v+=4){
	ldop_d->LdopDeflationBasisInit(v);
      }
      ldop_d->LdopDeflationBasisDiagonalise(LdopDeflVecs);
      
  }

  //close_comms();

}

void Dirac_BFM_HDCG_Wrapper::solve_HDCG(Fermion_t solution[2], Fermion_t source[2], double residual, int maxit, int cb){
  linop.residual = residual;
  linop.max_iter = maxit;
  
  Fermion_t src = linop.allocFermion();// this may be useless - use psi_h in the parent class
  Fermion_t Mtmp= linop.allocFermion();

  int GlobalSize[4];
  int SubGridSize[4];
  int NodeCoord[4];
  
  for (int mu =0; mu< 4; mu++){
    GlobalSize[mu]  = CommonPrms::instance()->global_size(mu);
    SubGridSize[mu] = CommonPrms::instance()->local_size(mu);
    NodeCoord[mu]   = Communicator::instance()->ipe(mu);
  }

  
  BFM_HDCG_Extend<float>  *ldop_f;
  ldop_f = new BFM_HDCG_Extend<float>(HDCGParams.Ls,HDCGParams.NumberSubspace,HDCGParams.Block,HDCGParams.Block,GlobalSize, SubGridSize, NodeCoord, &linop,&linop_single);
  ldop_f->CloneSubspace<double>(*ldop_d);
  ldop_f->SetParams(HDCGParams);
  
  ldop_d->PcgType = PcgAdef2;// not in the original wrapper
  ldop_f->PcgType = PcgAdef2f;

  int threads = linop.threads;

  // hack
  cb=Odd;//?????????????????? why have to force this? (works)



#pragma omp parallel
  {
#pragma omp for
    for(int i=0;i<threads;i++) {

      int me = linop.thread_barrier();

      // checkboardsource  - Using CGDiagonalMee preconditioning 
      linop.MooeeInv(source[1-cb],tmp,DaggerNo,1-cb);
      linop.Meo     (tmp,src,cb,DaggerNo);
      linop.axpy    (src,src,source[cb],-1.0);
      linop.MooeeInv(src,tmp,DaggerNo,cb);//my add CGdiagonalMee
      linop.Mprec   (tmp,src,Mtmp,DaggerYes);

      linop.ThreadBossLog("Solving for flexible = %d \n",HDCGParams.Flexible);
      ldop_f->AnalyseSpectrum=0; 

      if ( HDCGParams.Flexible == 1 ) { 
	linop_single.ThreadBossLog("Solving with single precision Ldop, fPcg\n");
	ldop_d->SloppyComms = 1; // DEBUG
	ldop_f->fPcg(ldop_d,solution[cb],src,tmp);
      } else if ( HDCGParams.Flexible == 2 ) {  
	linop.ThreadBossLog("Solving with double precision Ldop, fPcg\n");
	ldop_d->SloppyComms = 0;
	ldop_d->fPcg(ldop_d,solution[cb],src,tmp);
      } else { 
	linop_single.ThreadBossLog("Solving with single precision Ldop, Pcg\n");
	ldop_f->Pcg(ldop_d,solution[cb],src,tmp);
      }
      linop.thread_barrier();


      // DEBUG lines
      /*
      double sol_norm = sqrt(linop.norm(solution[cb]));
      double src_norm = sqrt(linop.norm(src));
      if (linop.isBoss() && !me ) printf("bfm_hdcg_wrapper: sol_norm   %10.8e \t cb = %d\n",sol_norm, cb);
      if (linop.isBoss() && !me ) printf("bfm_hdcg_wrapper: src_norm   %10.8e \t cb = %d\n",src_norm, cb);

      linop.Mprec(solution[cb],tmp,Mtmp,0);
      linop.Mprec(tmp,solution[1-cb],Mtmp,1); 
      linop.axpy(Mtmp,src,solution[1-cb],-1.0);
      double check_mprec = sqrt(linop.norm(Mtmp));
      if (linop.isBoss() && !me ) printf("bfm_hdcg_wrapper: recheck mprec    %10.8e \t cb = %d\n",check_mprec, cb);

      linop.Mprec(solution[cb],Mtmp,tmp,DaggerNo);// what is this doing????? no match with std mprecs MprecTilde???
      // Mprec test (comments assume cb=odd but valid for all cases)
      linop.Meo(solution[cb],tmp,1-cb, DaggerNo); // tmp   = Moe sol[cb]
      linop.MooeeInv(tmp,resid,DaggerNo,1-cb);   // resid = Mee-1 Moe sol[cb]
      linop.Meo(resid,tmp,cb, DaggerNo);          // tmp   = Meo Mee-1 Moe sol[cb]
      linop.MooeeInv(tmp,resid,DaggerNo,cb);     // resid = Moo-1 Meo Mee-1 Moe sol[cb]
      linop.axpy(tmp,resid,solution[cb],-1.0);   // tmp   = sol[cb] - Moo-1 Meo Mee-1 Moe sol[cb]
      
      linop.axpy(tmp,tmp,Mtmp,-1.0);               // tmp = Mtmp - tmp   should be zero

      double check_mprec1 = sqrt(linop.norm(tmp));
      if (linop.isBoss() && !me ) printf("bfm_hdcg_wrapper: mprec check      %10.8e \t cb = %d\n",check_mprec1, cb);
      */
      //////////////////////////////// end of debug lines



      // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
      linop.Meo(solution[cb],tmp,1-cb,DaggerNo);
      linop.axpy(src,tmp,source[1-cb],-1.0);
      linop.MooeeInv(src,solution[1-cb],DaggerNo,1-cb);

      // cb defect
      linop.Meo(solution[1-cb],tmp,cb,DaggerNo);
      linop.Mooee(solution[cb],Mtmp,DaggerNo,cb);
      linop.axpy(Mtmp,Mtmp,tmp,1.0);
      linop.axpy(Mtmp,Mtmp,source[cb],-1.0);
    
      double rf1 = linop.norm(Mtmp);
      //if (linop.isBoss() && !me ) printf("bfm_hdcg_wrapper: debug cb defect is %le\n",sqrt(rf1)); 


      // 1-cb defect
      linop.Meo(solution[cb],tmp,1-cb,DaggerNo);
      linop.Mooee(solution[1-cb],Mtmp,DaggerNo,1-cb);
      linop.axpy(Mtmp,Mtmp,tmp,1.0);
      linop.axpy(Mtmp,Mtmp,source[1-cb],-1.0);

      double rf=linop.norm(Mtmp);
      //if (linop.isBoss() && !me ) printf("bfm_hdcg_wrapper: debug 1-cb defect is %le\n",sqrt(rf)); // by construction = 0
      rf+= rf1;

      double f1 = linop.norm(source[cb]);
      double f2 = linop.norm(source[1-cb]);
      double f = f1+f2;

      if (linop.isBoss() && !me ) printf("bfm_hdcg_wrapper: unprec sol true residual is %le\n",sqrt(rf/f));

      /*
      //DEBUG lines
      if (linop.isBoss() && !me ) printf("bfm_hdcg_wrapper: debug source norm is %le\n",sqrt(f));

      double s1 = linop.norm(solution[cb]);
      double s2 = linop.norm(solution[1-cb]);

      if (linop.isBoss() && !me ) printf("bfm_hdcg_wrapper: debug source even norm is %le\n",sqrt(f2)); 
      if (linop.isBoss() && !me ) printf("bfm_hdcg_wrapper: debug source odd  norm is %le\n",sqrt(f1));
      if (linop.isBoss() && !me ) printf("bfm_hdcg_wrapper: debug sol    even norm is %le\n",sqrt(s2)); 
      if (linop.isBoss() && !me ) printf("bfm_hdcg_wrapper: debug sol    odd  norm is %le\n",sqrt(s1));
      */
 

    }
  }

  linop.freeFermion(src);
  linop.freeFermion(Mtmp);
  ldop_f->end();
  delete ldop_f;

}    


void Dirac_BFM_HDCG_Wrapper::solve_HDCG(FermionField &solution, FermionField &source){
  if (is_initialized) {
    // Wrapper of the previous basic method for the inverter
    int Nvol5d = Nvol_*BFMparams.Ls_;
    
    long double export_timing;
    FINE_TIMING_START(export_timing);
    BFM_interface.GaugeExport_to_BFM(u_);
    FINE_TIMING_END(export_timing);
    _Message(TIMING_VERB_LEVEL, "[Timing] - Dirac_BFM_Wrapper::solve_CGNE"
	     << " - Gauge Export to BFM timing = "
	     << export_timing << std::endl);

    // Initialize single precision gauge field if necessary
    if (HDCGParams.Flexible == 1 ){
      FINE_TIMING_START(export_timing);
      BFM_interface_single.GaugeExport_to_BFM(u_);
      FINE_TIMING_END(export_timing);
      _Message(TIMING_VERB_LEVEL, "[Timing] - Dirac_BFM_Wrapper::solve_CGNE"
	       << " - Gauge Export Single to BFM timing = "
	       << export_timing << std::endl);
      
      
    }

    int cb = Even;// by default
    
    // Load the fermion field to BFM
    // temporary source vector for the BFM data
    FermionField BFMsource(Nvol5d);
    FermionField BFMsol(Nvol5d);
    
    int half_vec = source.size()/(2*BFMparams.Ls_);
    for (int s =0 ; s< BFMparams.Ls_; s++){
      for (int i = 0; i < half_vec; i++){
	
	BFMsource.data.set(i+    2*s*half_vec, source.data[i+half_vec*s]);//even part
	BFMsource.data.set(i+(2*s+1)*half_vec, source.data[i+half_vec*(s+BFMparams.Ls_)]);//odd part 
	
	BFMsol.data.set(i+    2*s*half_vec, solution.data[i+half_vec*s]);//just even part is enough
	BFMsol.data.set(i+(2*s+1)*half_vec, solution.data[i+half_vec*(s+BFMparams.Ls_)]);//just even part is enough
      }
    }
    
    LoadSource(BFMsource,cb);
    LoadSource(BFMsource,1-cb);   
    LoadGuess(BFMsol, cb);
    LoadGuess(BFMsol, 1-cb);
    
    //////////////////////////// Execute the solver
    solve_HDCG(chi_h,psi_h, parameters.residual, parameters.max_iter, cb);
    ///////////////////////////////////////////////
    
    // Get the solution from BFM
    GetSolution(BFMsol,cb);
    GetSolution(BFMsol,1-cb);
    for (int s =0 ; s< BFMparams.Ls_; s++){
      for (int i = 0; i < half_vec; i++){
	solution.data.set(i+half_vec*s,                 BFMsol.data[i+2*s*half_vec]);//copy the even part
	solution.data.set(i+half_vec*(s+BFMparams.Ls_), BFMsol.data[i+(2*s+1)*half_vec]);//copy the odd part
      }
    }
  } else {
    std::ostringstream err_msg;
    err_msg << "Dirac_BFM_HDCG_Wrapper is not initialized";
    Errors::BaseErr("DEBUG", err_msg);
  }



}

Fermion_t* Dirac_BFM_HDCG_Wrapper::mult_inv_4d_base(Fermion_t psi_in[2]){
  // override the base class method with the same name
  if(is_initialized){
    BFM_interface.GaugeExport_to_BFM(u_);
    if (HDCGParams.Flexible == 1)
      BFM_interface_single.GaugeExport_to_BFM(u_);
    
    int cb = Even;

#pragma omp parallel for
    for (int t=0;t<threads;t++){
      linop.copy(psi_h[0], psi_in[0]);
      linop.copy(psi_h[1], psi_in[1]);
    }
    
    solve_HDCG(chi_h, psi_h, BFMparams.target_, BFMparams.max_iter_, cb);// cb is even 
    
  } else {
    CCIO::cout << "The operator was not initialized yet\n";
  }

}

