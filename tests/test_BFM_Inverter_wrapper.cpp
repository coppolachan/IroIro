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


using namespace std;

int Test_Solver_BFM::run(){
  CCIO::cout << ".::: Test Solver BFM Moebius" 
	     <<std::endl;


  Mapping::init_shiftField();

  //Parameters
  int Nvol  =  CommonPrms::instance()->Nvol();
  double mq = 0.0;
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

  Dirac_optimalDomainWall_EvenOdd DWF_EO(b,c, M5, mq, omega, &Dw_eo, &Dw_oe, &(conf_.data));
  
  ////////////////////////////////////////////////////////////
  //    Setup DWF operator
  ////////////////////////////////////////////////////////////
  Dirac_BFM_Wrapper BFMoperator(&conf_.data);
  BFMoperator.set_ScaledShamirCayleyTanh(mq, M5, Ls, ht_scale);
  BFMoperator.solver_params(50000, 1e-12);
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
  CCIO::cout << "phi,size,norm: " << phi.size() << " " << phi_norm << endl;

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
  ///////////////////////////////////////////////////////////////////////////////////////
  // Solver using internal Dirac_optimalDomainWall_EvenOdd method solve_eo
  SolverOutput SO;
  Field output_f(vphi);
  DWF_EO.solve_eo(output_f,fe, SO,  10000, 1.0e-24);
  SO.print();
  
  ///////////////////////////////////////////////////////////////////////////////////////
 // Solver using BFM
  FermionField BFMsolution(Nvol5d);
  BFMoperator.solve_CGNE(BFMsolution, EO_source);
  
  int vect4d_hsize = fe.size()/Ls;   
  for (int s =0 ; s< Ls; s++){
    for (int i = 0; i < vect4d_hsize; i++){
      fe.set(i+vect4d_hsize*s, BFMsolution.data[i+2*s*vect4d_hsize]);
    }
  }
  //////////////////////////////////////////////////////////////////////////////////  
  // Check
  Field F_Diff;
  F_Diff = output_f;
  F_Diff -= fe;
  /*
  for (int i = 0; i < fe.size() ; i++){
    double diff = abs(fe[i]-output_f[i]);
    if (diff>1e-8) CCIO::cout << "*";
    CCIO::cout << "["<<i<<"] "<<fe[i] << "  "<<output_f[i]
               << "  "<< diff << "\n";
  }
  */

  CCIO::cout << "Operator Difference BFM-IroIro = "<< F_Diff.norm() << "\n";
  
  return 0;
}
