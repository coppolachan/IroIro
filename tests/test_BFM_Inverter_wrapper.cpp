/*!
 * @file test_Solver_BFM.cpp
 * @brief Definition of classes for testing the BFM classes
 */
#include "test_Solver_BFM.hpp"
#include "Dirac_ops/dirac_wilson_EvenOdd.hpp"
#include "Geometry/siteIndex_EvenOdd.hpp"
#include "Measurements/FermionicM/quark_prop_meas_factory.hpp"
#include "Measurements/FermionicM/qprop_DomainWall.hpp"
#include "Measurements/FermionicM/source_types.hpp"

#include "Dirac_ops/BFM_Wrapper/dirac_BFM_wrapper.hpp"
#include "Dirac_ops/BFM_Wrapper/dirac_BFM_wrapper_factory.hpp"

#include "Action/action_fermiontype_factory.hpp"
#include "Tools/Architecture_Optimized/utils_BGQ.hpp"


using namespace std;

int Test_Solver_BFM::run(){
  CCIO::cout << ".::: Test Solver BFM Factories" 
	     <<std::endl;


  Mapping::init_shiftField();

  //Parameters
  int Nvol  =  CommonPrms::instance()->Nvol();
  double mq = 0.01;
  double M5 = -1.6;
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
  Dirac_Wilson Dw_eo(M5,&(conf_.data),Dop::EOtag());
  Dirac_Wilson Dw_oe(M5,&(conf_.data),Dop::OEtag());

  InputConfig Iconf(&conf_);

  Dirac_DomainWall_EvenOdd DWF_EO(b,c, M5, mq,
					 omega, &Dw_eo, &Dw_oe);
  
  ////////////////////////////////////////////////////////////
  //    Setup DWF operator
  ////////////////////////////////////////////////////////////
  XML::node top_Node = DWFnode;
  XML::node SolverNode = DWFnode;
  XML::descend (SolverNode, "Solver",MANDATORY);

  CCIO::cout << "Creating factory\n";
  DiracBFMoperatorFactory BFMfactory(DWFnode);
  CCIO::cout << "Creating object\n";
  
  Dirac_BFM_Wrapper* BFMop_with_fact = BFMfactory.getDirac(Iconf);
  BFMop_with_fact->set_SolverParams(SolverNode);
  BFMop_with_fact->initialize();
  
  Dirac_BFM_Wrapper BFMoperator(DWFnode,&conf_.data, &DWF_EO);// standard creation
  BFMoperator.set_SolverParams(SolverNode);
  BFMoperator.initialize();

  //////////////////////////////////////////////////////////// 
  // Sources setup
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
  CCIO::cout << "fe norm : " << fe.norm() << endl;

  FermionField EO_source(Nvol5d);

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
  FermionField fe_FF(fe); 
  ////////////////////////////////////////////////////////////////////////////////// 
  // Force term 
  Field force_BFM = BFMoperator.md_force(fe ,fe);
  Field force_IroIro =  DWF_EO.md_force(fe, fe);


  Field ForceDiff =  force_BFM;
  ForceDiff -= force_IroIro;
  CCIO::cout << "Force check BFM - IroIro  "<< force_BFM.norm() 
  	     << "  "<< force_IroIro.norm() << "\n";

  CCIO::cout << "Force check Difference BFM-IroIro = "<< ForceDiff.norm() << "\n";

  ///////////////////////////////////////////////////////////////////////////////////////
  // Solver CG using SolverFactory
  SolverCG_DWF_opt_Factory SolveBFM(SolverNode); 
  Solver_CG_DWF_Optimized *SolverCGNE = SolveBFM.getSolver(BFMop_with_fact);
  FermionField BFMsolution3(Nvol5d/2);//half vector
  SolverCGNE->solve(BFMsolution3.data,fe);


  // Using the factory object
  FermionField BFMsolution2(Nvol5d);
  BFMop_with_fact->solve_CGNE(BFMsolution2, fe_FF);
  
  // Using the factory object (mixed precision)
  FermionField BFMsolutionMixed(Nvol5d);
  BFMop_with_fact->solve_CGNE_mixed_prec(BFMsolutionMixed, fe_FF);

  ///////////////////////////////////////////////////////////////////////////////////////
  // Action Nf2 using BFM
  SmartConf SConf;
  SConf.ThinLinks = &conf_;
  
  XML::descend(top_Node, "TestAction", MANDATORY);
  TwoFlavorDomainWall5dEO_BFM_ActionFactory Act(top_Node);
  Act.getAction(&conf_, &SConf);
  

  ///////////////////////////////////////////////////////////////////////////////////////
  // Solver using internal Dirac_DomainWall_EvenOdd method solve_eo
  SolverOutput SO;
  Field output_f(vphi);
  DWF_EO.solve_eo(output_f,fe, SO,  10000, 1.0e-24);
  SO.print();

  ///////////////////////////////////////////////////////////////////////////////////////
  // Solver MultishiftCG using BFM
  vector_double shifts(5);
  shifts[0] = 0.001;
  shifts[1] = 0.002;
  shifts[2] = 0.003;
  shifts[3] = 0.004;
  shifts[4] = 0.005;
  vector_double mresiduals(5);
  mresiduals[0] = 1e-12;
  mresiduals[1] = 1e-12;
  mresiduals[2] = 1e-12;
  mresiduals[3] = 1e-12;
  mresiduals[4] = 1e-12;
  
  vector_Field shifted_sol;
  shifted_sol.resize(shifts.size());
  for (int i=0; i< shifted_sol.size(); ++i) {
    shifted_sol[i].resize(fe.size());
  }

  DWF_EO.solve_ms_eo(shifted_sol,fe, SO, shifts,  10000, 1.0e-24);
  SO.print();

  std::vector < FermionField > BFM_ms_solution(shifts.size());
 
  
  BFMop_with_fact->solve_CGNE_multishift(BFM_ms_solution, fe_FF, shifts, mresiduals);
  /*
  for(int nsol = 0; nsol < 5; nsol++){
    Field F_Diff = shifted_sol[nsol];
    //F_Diff -= fe;
    F_Diff -= BFM_ms_solution[nsol].data;
    CCIO::cout << "Operator IroIro ["<<nsol<<"] = "<< shifted_sol[nsol].norm() << "\n";
    CCIO::cout << "Operator BFM ["<<nsol<<"] = "<< BFM_ms_solution[nsol].norm() << "\n";
 
    CCIO::cout << "Operator Difference BFM-IroIro ["<<nsol<<"] = "<< F_Diff.norm() << "\n";
  }
  */
  ///////////////////////////////////////////////////////////////////////////////////////
  // Solver MultishiftCG using BFM in mixed precision
  BFMop_with_fact->solve_CGNE_multishift_mixed_precision(BFM_ms_solution, fe_FF, shifts, mresiduals);

 
  for(int nsol = 0; nsol < 5; nsol++){
    Field F_Diff = shifted_sol[nsol];
    //F_Diff -= fe;
    F_Diff -= BFM_ms_solution[nsol].data;
    CCIO::cout << "Operator IroIro ["<<nsol<<"] = "<< shifted_sol[nsol].norm() << "\n";
    CCIO::cout << "Operator BFM ["<<nsol<<"] = "<< BFM_ms_solution[nsol].norm() << "\n";
 
    CCIO::cout << "Operator Difference BFM-IroIro ["<<nsol<<"] = "<< F_Diff.norm() << "\n";
  }
  

  /*
  Field sol(fe);
  Spinor *fe_ptr = (Spinor*)fe.getaddr(0);
  Spinor *sol_ptr = (Spinor*)sol.getaddr(0);
  
  
  //sol = source;
  //sol *= ConstTerm; 
  for (int i = 0; i < Ls; i++)
    BGWilsonLA_MultScalar(sol_ptr+i*CommonPrms::Nvol(), 
			  fe_ptr+i*CommonPrms::Nvol(), 
			  2.0, CommonPrms::Nvol());
  


  ////////////////////////////////////////////////////////////////////////////////// 
  // Mult & Mult_dag
  
  Field mult_BFM = BFMoperator.mult(fe);
  Field mult_IroIro = DWF_EO.mult(fe);

  Field multDiff =  mult_BFM;
  multDiff -= mult_IroIro;
  CCIO::cout << "Mult check BFM - IroIro  "<< mult_BFM.norm() 
  	     << "  "<< mult_IroIro.norm() << "\n";

  CCIO::cout << "Mult check Difference BFM-IroIro = "<< multDiff.norm() << "\n";

  Field mult_dag_BFM = BFMoperator.mult_dag(fe);
  Field mult_dag_IroIro = DWF_EO.mult_dag(fe);

  Field mult_dag_Diff =  mult_dag_BFM;
  mult_dag_Diff -= mult_dag_IroIro;
  CCIO::cout << "Mult_dag check BFM - IroIro  "<< mult_dag_BFM.norm() 
  	     << "  "<< mult_dag_IroIro.norm() << "\n";

  CCIO::cout << "Mult_dag check Difference BFM-IroIro = "<< mult_dag_Diff.norm() << "\n";
  /*
  
  ///////////////////////////////////////////////////////////////////////////////////////
  // Solver CG using BFM
  
  FermionField BFMsolution(Nvol5d);
  BFMoperator.solve_CGNE(BFMsolution, EO_source);
  */

  
 
  Field F_Diff = output_f;
  //F_Diff -= BFMsolution3.data;
  //F_Diff -= BFMsolution2.data;
  F_Diff -= BFMsolutionMixed.data;
  /*
  for (int i = 0; i < fe.size() ; i++){
    double diff = abs(BFMsolution3.data[i]-output_f[i]);
    if (diff>1e-8) CCIO::cout << "*";
    CCIO::cout << "["<<i<<"] "<<BFMsolution3.data[i] << "  "<<output_f[i]
               << "  "<< diff << "\n";
  }
  */
  CCIO::cout << "Operator Difference BFM-IroIro = "<< F_Diff.norm() << "\n";
  
  return 0;
}
