/*!
 * @file test_Solver_HDCH.cpp
 * @brief Definition of classes for testing the BFM classes
 */
#include "test_Solver_HDCG.hpp"
#include "Dirac_ops/dirac_wilson_EvenOdd.hpp"
#include "Geometry/siteIndex_EvenOdd.hpp"
#include "Dirac_ops/BFM_Wrapper/dirac_BFM_wrapper_factory.hpp"
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

  InputConfig Iconf(&conf_); 
  Mapping::init_shiftField();

  //Parameters
  int Nvol  =  CommonPrms::instance()->Nvol();
  double mq = 0.001;
  double M5 = 1.6;
  int Ls    = 8;
  double ht_scale = 2.0;
  int Nvol5d = Nvol*Ls;
  double c   = 1.0;
  double b   = 2.0;
  vector<double> omega(Ls,1.0);

  double DWFresidual = 1e-24; //square of bfm residual
  double DWFmax_iter = 2000;


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

  // creation of Dirac_Wilson kernel operators 
  Dirac_Wilson Dw_eo(-M5,&(conf_.data),Dop::EOtag());
  Dirac_Wilson Dw_oe(-M5,&(conf_.data),Dop::OEtag());
  Dirac_DomainWall_EvenOdd DWF_EO(b,c, -M5, mq, omega, &Dw_eo, &Dw_oe);
  
  /********************************************************

   * Setup HDCG operator

   ********************************************************
   */

  // Smoother control?
    
  //////////////////////////////////////
  // Create and initializes the HDCG class
  XML::node SolverNode = DWFnode;
  XML::descend (SolverNode, "Solver",MANDATORY);

  /*
  //Standard creation
  Dirac_BFM_HDCG_Wrapper HDCG_Solver(DWFnode,       // XML node describing the operator
				     &conf_.data,   // configuration
				     &DWF_EO);      // passes the IroIro operator
  */

  CCIO::cout << "Creating BFM factory\n";
  DiracBFMoperatorFactory BFMfactory(DWFnode);
  CCIO::cout << "Creating BFM object\n";

  Dirac_BFM_Wrapper* BFMop_with_fact = BFMfactory.getDirac(Iconf);
  BFMop_with_fact->set_SolverParams(SolverNode);
  BFMop_with_fact->initialize();


  // Test factory creation
  CCIO::cout << "Creating HDCG factory\n";
  DiracBFM_HDCGoperatorFactory BFM_HDCGfactory(DWFnode);
  CCIO::cout << "Creating HDCG object\n";

  Dirac_BFM_HDCG_Wrapper* HDCG_Solver_with_fact = BFM_HDCGfactory.getDirac(Iconf);
  CCIO::cout << "Created HDCG object\n";

  HDCG_Solver_with_fact->HDCG_init(DWFnode);  // Initialization (in dirac_BFM_HDCG.cpp)
  CCIO::cout << "Tester Msg: HDCG_Solver.HDCG_init completed \n";

  HDCG_Solver_with_fact->set_SolverParams(SolverNode);
  CCIO::cout << "Tester Msg: HDCG_Solver.set_SolverParams completed \n";
  
  //////////////////////////////////////////////
  // launch the test
  // First set the source vectors
  FermionField source(src.mksrc(0,0));
 
  SiteIndex_EvenOdd* ieo = SiteIndex_EvenOdd::instance();
  vector<int> esec= ieo->esec();
  vector<int> osec= ieo->osec();
  // source generation random 

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
  
  // copy fields in a single eo source
  assert(fe.size() == fo.size());
  int vect4d_hsize=fe.size()/Ls;
  for (int s =0 ; s< Ls; s++){
    for (int i = 0; i < fe.size()/Ls; i++){
      EO_source.data.set(i+vect4d_hsize*s,      fe[i+vect4d_hsize*s]);
      EO_source.data.set(i+vect4d_hsize*(s+Ls), fo[i+vect4d_hsize*s]);
    }
  }
  double EO_source_norm = EO_source.norm();
  CCIO::cout << "EO_source filled - norm:  "<< EO_source_norm << "\n";
  CCIO::cout << "EO_source filled - Even norm:  "<< fe.norm() << "\n";
  CCIO::cout << "EO_source filled - Odd  norm:  "<< fo.norm() << "\n";
  FermionField fe_FF(fe);

  ///////////////////////////////////////////////////////////////////////////////////////
  FermionField BFMSolution(Nvol5d);
  FermionField BFMSolutionMixed(Nvol5d);

  // Test solver
  BFMop_with_fact->solve_CGNE_mixed_prec(BFMSolutionMixed,fe_FF);


  // 1. Init subspace (computes the vectors)
  CCIO::cout << ".::::::::::::::  Init subspace\n";
  HDCG_Solver_with_fact->HDCG_subspace_init();
  // 2. Compute the subspace (coarse space Little Dirac operator construction)
  CCIO::cout << ".::::::::::::::  Compute subspace\n";
  HDCG_Solver_with_fact->HDCG_subspace_compute(0);  // 0 = double, 1 = single

  // 3. Refine the subspace (optional)
  if (HDCG_Solver_with_fact->do_refine()){
    CCIO::cout << ".::::::::::::::  Refine subspace\n";
    HDCG_Solver_with_fact->HDCG_subspace_refine(); 
    HDCG_Solver_with_fact->HDCG_subspace_compute(0);  // 0 = double, 1 = single
  }

  // 4. Actually solve with CG (solves M x = b , not MdagM x = b)
  CCIO::cout << ".::::::::::::::  Solve using CG\n";
  HDCG_Solver_with_fact->solve_HDCG(BFMSolution,EO_source);
  
  // 5. Free pointers
  CCIO::cout << "Subspace free pointers\n";
  HDCG_Solver_with_fact->HDCG_subspace_free();
  ////////////////////////////////////// end of HDCG

  // Solver using internal Dirac_DomainWall_EvenOdd method solve_eo
  SolverOutput SO;
  Field output_f(vphi);
  // Must precondition the source because the IroIro solver is giving the solution of
  // MdagM chi = psi
  //
  // Preconditioned source in BFM terminology
  // prec_fe = Mprec^dag (Mee^-1(fe - MoeMoo^-1 fo) 

  Field prec_tmp1 = DWF_EO.mult_ee_inv(fe);
  prec_tmp1 -= DWF_EO.mult_eo(DWF_EO.mult_oo_inv(fo));
  Field prec_fe = DWF_EO.mult_dag(prec_tmp1);

  // Solve with IroIro (reference result) - Solving MdagM x = b 
  DWF_EO.solve_eo(output_f,prec_fe, SO,  DWFmax_iter, DWFresidual);
  SO.print();
  
  // Check IroIro preconditioned solution even - MdagM
  Field IroIroSolveCheck = (DWF_EO.mult(output_f));
  IroIroSolveCheck -= prec_tmp1;
  CCIO::cout << "IroIro MdagM check = "<< IroIroSolveCheck.norm() << "\n";

  // Get the Odd part of the unpreconditioned solution
  prec_tmp1 = fo;
  prec_tmp1 -= DWF_EO.mult_oo(DWF_EO.mult_oe(output_f));
  Field IroIroOdd = DWF_EO.mult_oo_inv(prec_tmp1);
  
  //Test1  :  Mee sol_e + Moe sol_o - src_e = 0 
  Field test1 = DWF_EO.mult_ee(output_f);
  test1 += DWF_EO.mult_ee(DWF_EO.mult_eo(IroIroOdd));
  test1 -= fe;
  CCIO::cout << "IroIro unpreconditioned M check Even  = "<< test1.norm() << "\n";

  //Test2  :  Moo sol_o + Meo sol_e - src_o = 0 
  test1 = DWF_EO.mult_oo(IroIroOdd);
  test1 += DWF_EO.mult_oo(DWF_EO.mult_oe(output_f));
  test1 -= fo;
  CCIO::cout << "IroIro unpreconditioned M check Odd   = "<< test1.norm() << "\n";


  //////////////////////////////////////////////////////////////////////////////////
  // Separate the even and odd parts of the BFM solutions to different fields
  Field bfm_se(fe);
  Field bfm_so(fo);
  for (int s =0 ; s< Ls; s++){
    for (int i = 0; i < vect4d_hsize; i++){
      bfm_se.set(i+vect4d_hsize*s, BFMSolution.data[i+(s   )*vect4d_hsize]);
      bfm_so.set(i+vect4d_hsize*s, BFMSolution.data[i+(s+Ls)*vect4d_hsize]);
    }
  }

  CCIO::cout << "BFM sol even norm = "<< bfm_se.norm() << "\n";
  CCIO::cout << "BFM sol odd  norm = "<< bfm_so.norm() << "\n";

  // Check BFM solution against IroIro (preconditioned)
  Field Difference = output_f;
  Difference -= bfm_se;
  CCIO::cout << "Operator Difference BFM-IroIro (Even) = "<< Difference.norm() << "\n";

  Difference = IroIroOdd;
  Difference -= bfm_so;
  CCIO::cout << "Operator Difference BFM-IroIro (Odd)  = "<< Difference.norm() << "\n";


  

  
  return 0;
}
