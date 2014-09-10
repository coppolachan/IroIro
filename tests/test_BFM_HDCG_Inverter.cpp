/*!
 * @file test_Solver_HDCH.cpp
 * @brief Definition of classes for testing the BFM classes
 */
#include "test_Solver_HDCG.hpp"
#include "Dirac_ops/dirac_wilson_EvenOdd.hpp"
#include "Geometry/siteIndex_EvenOdd.hpp"
#include "Measurements/FermionicM/quark_prop_meas_factory.hpp"
#include "Measurements/FermionicM/qprop_DomainWall.hpp"
#include "Measurements/FermionicM/source_types.hpp"

//wrap
#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION
#include "bfm.h"
#include "BfmHDCG.h"
#include "Tools/Bagel/bfm_storage.hpp"
#include "Dirac_ops/BFM_Wrapper/dirac_BFM_HDCG.hpp"

using namespace std;

int Test_Solver_HDCG::run(){
  CCIO::cout << ".::: Test Solver HDCG BFM Moebius" 
	     <<std::endl;

 
  Mapping::init_shiftField();

  //Parameters
  int Nvol  =  CommonPrms::instance()->Nvol();
  double mq = 0.0;
  double M5 = 1.6;
  int Ls    = 8;
  double ht_scale = 2.0;
  int Nvol5d = Nvol*Ls;
  double c   = 1.0;
  double b   = 2.0;
  vector<double> omega(Ls,1.0);

  // Random number generator initialization
  unsigned long init[4]={0x123, 0x234, 0x345, 0x456};
  int length = 4;
  RandNum_MT19937 rand(init, length);
  
  vector<int> spos(4,0);
  spos[0] = 0;
  spos[1] = 0;
  spos[2] = 0;
  spos[3] = 0;

  double alpha = 5.0;
  Source_local<Format::Format_F> src(spos,Nvol5d);
  // Source_Gauss<Format::Format_F> src(spos,alpha, Nvol5d);

  // Testing the correctness of the operator

  // creation of Dirac_Wilson kernel operators 
  Dirac_Wilson Dw_eo(-M5,&(conf_.data),Dop::EOtag());
  Dirac_Wilson Dw_oe(-M5,&(conf_.data),Dop::OEtag());


  Dirac_DomainWall_EvenOdd DWF_EO(b,c, -M5, mq, omega, &Dw_eo, &Dw_oe);
  
  /********************************************************

   * Setup DWF operator

   ********************************************************
   */

  bfmarg  dwfa;
  dwfa.node_latt[0]  = CommonPrms::instance()->Nx();
  dwfa.node_latt[1]  = CommonPrms::instance()->Ny();
  dwfa.node_latt[2]  = CommonPrms::instance()->Nz();
  dwfa.node_latt[3]  = CommonPrms::instance()->Nt();

  for(int mu=0;mu<4;mu++){
    dwfa.neighbour_plus[mu]  = Communicator::instance()->node_up(mu);
    dwfa.neighbour_minus[mu] = Communicator::instance()->node_dn(mu);
  }
  dwfa.verbose = 8;
  dwfa.time_report_iter=100;

  for(int mu=0;mu<4;mu++){
    if ( (CommonPrms::instance()->NPE(mu))>1 ) {
      dwfa.local_comm[mu] = 0;
    } else {
      dwfa.local_comm[mu] = 1;
    }
  }

  int threads = 64;
  bfmarg::Threads(threads);

  bfmarg::UseCGdiagonalMee(1);


  dwfa.ScaledShamirCayleyTanh(mq,M5,Ls,ht_scale);
  dwfa.rb_precondition_cb=Even;
  dwfa.max_iter=50000;
  dwfa.residual=1.0e-12;


  // Define the HDCG parameters 
  BfmHDCGParams HDCGParams;

  //  HDCGcontrol Control = HdcgGenerateSubspace;////?
  HDCGParams.NumberSubspace = 24;
  HDCGParams.Ls = Ls;
  HDCGParams.Block[0] = 4;
  HDCGParams.Block[1] = 4;
  HDCGParams.Block[2] = 4;
  HDCGParams.Block[3] = 4;
  HDCGParams.Block[4] = Ls;

  HDCGParams.SubspaceRationalRefine = 0;
  HDCGParams.SubspaceRationalRefineLo =   0.003;
  HDCGParams.SubspaceRationalRefineResidual =   1.0e-3;
  HDCGParams.SubspaceRationalLs = 6;
  HDCGParams.SubspaceRationalLo =   0.001;
  HDCGParams.SubspaceRationalMass =   0.002;
  HDCGParams.SubspaceRationalResidual =   1.0e-4;
  HDCGParams.SubspaceSurfaceDepth = 6;
  HDCGParams.LittleDopSolverResidualInner =   0.04;
  HDCGParams.LittleDopSolverResidualVstart =   0.04;
  HDCGParams.LittleDopSolverResidualSubspace =   1.0e-10;
  HDCGParams.LittleDopSubspaceRational = 0;
  HDCGParams.LittleDopSolverIterMax = 10000;
  HDCGParams.LdopDeflVecs = 24;
  HDCGParams.PreconditionerKrylovResidual =   1.0e-5;
  HDCGParams.PreconditionerKrylovIterMax = 7;
  HDCGParams.PreconditionerKrylovShift =   1.0;
  HDCGParams.PcgSingleShift =   0.0;
  HDCGParams.LittleDopSolver = 0;
  HDCGParams.Flexible = 1;//?
 


  //////////////////////////////////////
  // Create the HDCG class
  XML::node SolverNode = DWFnode;
  XML::descend (SolverNode, "Solver",MANDATORY);

  Dirac_BFM_HDCG_Wrapper HDCG_Solver(DWFnode,&conf_.data, &DWF_EO);
  HDCG_Solver.HDCG_init(HDCGParams, dwfa);
  HDCG_Solver.set_SolverParams(SolverNode);
  HDCG_Solver.initialize(); // needs the operator and solver params to be set

  
  //launch the test
  FermionField source(src.mksrc(0,0));
 
  SiteIndex_EvenOdd* ieo = SiteIndex_EvenOdd::instance();
  vector<int> esec= ieo->esec();
  vector<int> osec= ieo->osec();
  // source generation (spin, color)
  // Field fe = src.mksrc(esec,2,0);
  // Field fo = src.mksrc(osec,2,0);
  valarray<double> vphi(DWF_EO.fsize());
  rand.get(vphi);
  Field phi(vphi);    // phi: generated from random numbers
  rand.get(vphi);
  Field fe(vphi);
  rand.get(vphi);
  Field fo(vphi);
  double phi_norm = phi.norm();
  CCIO::cout << "phi,size,norm: " << phi.size() << " " << phi_norm << endl;

#ifdef DEBUG
  valarray<double> phi(source.format.size()/2);
  MPrand::mp_get_gauss(phi, rand);
  Field fe(phi);
  Field fo(phi);
#endif
  FermionField EO_source(Ls*Nvol);
  CCIO::cout << "EO_source size: "<< EO_source.size() <<"\n";
  
  // copy fields in the eo source CB blocked
  // interleaving the e/o blocks
  assert(fe.size() == fo.size());
  for (int s =0 ; s< Ls; s++){
    for (int i = 0; i < fe.size()/Ls; i++){
      EO_source.data.set(i+2*fe.size()/Ls*s, fe[i+fe.size()/Ls*s]);
      EO_source.data.set(i+(2*s+1)*fe.size()/Ls, fo[i+fe.size()/Ls*s]);
    }
  }
  CCIO::cout << "EO_source filled\n";
  ///////////////////////////////////////////////////////////////////////////////////////
  FermionField BFMSolution(Nvol5d);
  // Test CGNE solver

  // 1. Init subspace 
  CCIO::cout << "Init subspace\n";
  HDCG_Solver.HDCG_subspace_init();
  // 2. Compute the subspace (coarse space Little Dirac operator construction)
  CCIO::cout << "Compute subspace\n";
  HDCG_Solver.HDCG_subspace_compute(0);
    
  // 3. Actually solve with CG
  CCIO::cout << "Solve using CG\n";
  HDCG_Solver.solve_HDCG(BFMSolution,EO_source, 1.0e-12, 50000);
  
  // 4. Free pointers
  CCIO::cout << "Subspace free pointers\n";
  HDCG_Solver.HDCG_subspace_free();


  // Solver using internal Dirac_DomainWall_EvenOdd method solve_eo
  SolverOutput SO;
  Field output_f(vphi);
  DWF_EO.solve_eo(output_f,fe, SO,  10000, dwfa.residual*dwfa.residual);
  SO.print();
  
  //Apply operator back on the inverse
  FermionField IroIroFull(Nvol5d);

  Field IroIroSol_even = DWF_EO.mult_dag(DWF_EO.mult(fe));
  Field IroIroSol_odd  = DWF_EO.mult_dag(DWF_EO.mult(fo));

  /*
  for (int i = 0; i < IroIroFull.size() ; i++){
    double diff = abs(IroIroFull.data[i]-EO_source.data[i]);
    if (diff>1e-8) CCIO::cout << "*";
    CCIO::cout << "["<<i<<"] "<<IroIroFull.data[i] << "  "<<EO_source.data[i]
               << "  "<< diff << "\n";
  }
  */
  
  FermionField Difference(Nvol5d);
  Difference = EO_source;
  Difference -= BFMSolution;

  /*
  Field F_Diff;
  F_Diff = output_f;
  F_Diff -= fe;
  
  for (int i = 0; i < fe.size() ; i++){
    double diff = abs(fe[i]-output_f[i]);
    if (diff>1e-8) CCIO::cout << "*";
    CCIO::cout << "["<<i<<"] "<<fe[i] << "  "<<output_f[i]
               << "  "<< diff << "\n";
  }
  */
  CCIO::cout << "Operator Difference BFM-IroIro = "<< Difference.norm() << "\n";
  
  return 0;
}
