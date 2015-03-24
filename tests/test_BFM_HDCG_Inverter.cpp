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
#include <omp.h>


using namespace std;

int Test_Solver_HDCG::run(){
  CCIO::cout << ".::: Test Solver HDCG BFM Moebius" 
	     <<std::endl;

 
  Mapping::init_shiftField();

  //Parameters
  int Nvol  =  CommonPrms::instance()->Nvol();
  double mq = 0.01;
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
  dwfa.time_report_iter = 100;

  for(int mu=0;mu<4;mu++){
    if ( (CommonPrms::instance()->NPE(mu))>1 ) {
      dwfa.local_comm[mu] = 0;
    } else {
      dwfa.local_comm[mu] = 1;
    }
  }

  int threads = omp_get_max_threads();
  bfmarg::Threads(threads);
  bfmarg::UseCGdiagonalMee(1);

  dwfa.ScaledShamirCayleyTanh(mq,M5,Ls,ht_scale);
  dwfa.rb_precondition_cb = Even;
  dwfa.max_iter           = 50000;
  dwfa.residual           = 1.0e-12;


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

  HDCGParams.SubspaceRationalRefine = 0;// single /double? 
  HDCGParams.SubspaceRationalRefineLo =   0.02;
  HDCGParams.SubspaceRationalRefineResidual =   1.0e-3;
  HDCGParams.SubspaceRationalLs = 8;
  HDCGParams.SubspaceRationalLo =   0.02;
  HDCGParams.SubspaceRationalMass =   0.001;
  HDCGParams.SubspaceRationalResidual =   1.0e-6;
  HDCGParams.SubspaceSurfaceDepth = HDCGParams.SubspaceRationalLs;
  HDCGParams.LittleDopSolverResidualInner =   0.04;
  HDCGParams.LittleDopSolverResidualVstart =   0.04;
  HDCGParams.LittleDopSolverResidualSubspace =   1.0e-6;
  HDCGParams.LittleDopSubspaceRational = 0;
  HDCGParams.LittleDopSolverIterMax = 100;
  HDCGParams.LdopDeflVecs = 24;
  HDCGParams.PreconditionerKrylovResidual =   1.0e-5;
  HDCGParams.PreconditionerKrylovIterMax = 7;
  HDCGParams.PreconditionerKrylovShift =   1.0;
  HDCGParams.PcgSingleShift =   0.0;
  HDCGParams.LittleDopSolver = LittleDopSolverCG;
  HDCGParams.Flexible = 2;// 2= double inner,  double outer
 


  //////////////////////////////////////
  // Create and initializes the HDCG class
  XML::node SolverNode = DWFnode;
  XML::descend (SolverNode, "Solver",MANDATORY);

  Dirac_BFM_HDCG_Wrapper HDCG_Solver(DWFnode,       // XML node describing the operator
				     &conf_.data,   // configuration
				     &DWF_EO);      // passes the IroIro operator
  HDCG_Solver.HDCG_init(HDCGParams, dwfa);  // Initialization (in dirac_BFM_HDCG.cpp)
  CCIO::cout << "HDCG_Solver.HDCG_init completed \n";
  

  HDCG_Solver.set_SolverParams(SolverNode);


  //HDCG_Solver.initialize(); // needs the operator and solver params to be set

  //////////////////////////////////////////////
  //launch the test
  // First set the source vectors
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
  double EO_source_norm = EO_source.norm();
  CCIO::cout << "EO_source filled - norm:  "<< EO_source_norm << "\n";

  ///////////////////////////////////////////////////////////////////////////////////////
  FermionField BFMSolution(Nvol5d);
  // Test solver

  // 1. Init subspace 
  CCIO::cout << ".::::::::::::::  Init subspace\n";
  HDCG_Solver.HDCG_subspace_init();
  // 2. Compute the subspace (coarse space Little Dirac operator construction)
  CCIO::cout << ".::::::::::::::  Compute subspace\n";
  HDCG_Solver.HDCG_subspace_compute(0);
    
  // 3. Actually solve with CG
  CCIO::cout << ".::::::::::::::  Solve using CG\n";
  HDCG_Solver.solve_HDCG(BFMSolution,EO_source, dwfa.residual , dwfa.max_iter);
  
  // 4. Free pointers
  CCIO::cout << "Subspace free pointers\n";
  HDCG_Solver.HDCG_subspace_free();


  // Solver using internal Dirac_DomainWall_EvenOdd method solve_eo
  SolverOutput SO;
  Field output_f(vphi);
  //precondition fe
  fe = (DWF_EO.mult_dag(fe));
  DWF_EO.solve_eo(output_f,fe, SO,  dwfa.max_iter, dwfa.residual*dwfa.residual);
  SO.print();
  

  // Copy the even and odd parts of the BFM solutions to different fields
  int vect4d_hsize=fe.size()/Ls;
  for (int s =0 ; s< Ls; s++){
    for (int i = 0; i < vect4d_hsize; i++){
      fe.set(i+vect4d_hsize*s, BFMSolution.data[i+(s   )*vect4d_hsize]);
      fo.set(i+vect4d_hsize*s, BFMSolution.data[i+(s+Ls)*vect4d_hsize]);
    }
  }

  // Check BFM solution
  Field IroIroSol_even = (DWF_EO.mult(fe));
  Field IroIroSol_odd  = (DWF_EO.mult(fo));


  Field Difference = output_f;
  Difference -= fe;

  /*
  for (int i = 0; i < fe.size()/200 ; i++){
    double diff = abs(fe[i]-output_f[i]);
    if (diff>1e-8) CCIO::cout << "*";
    CCIO::cout << "["<<i<<"] "<<fe[i] << "  "<<output_f[i]
               << "  "<< diff << "\n";
  }
  */

  /*
  for (int i = 0; i < fe.size()/200 ; i++){
    double diff = abs(EO_source.data[i]-IroIroSol_odd[i]);
    if (diff>1e-8) CCIO::cout << "*";
    CCIO::cout << "["<<i<<"] "<<EO_source.data[i] << "  "<<IroIroSol_odd[i]
               << "  "<< diff << "\n";
  }
  */

  CCIO::cout << "Operator Difference BFM-IroIro (Even) = "<< Difference.norm() << "\n";


  Difference = output_f;
  Difference -= fo;
  CCIO::cout << "Operator Difference BFM-IroIro (Odd)  = "<< Difference.norm() << "\n";

  
  return 0;
}
