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

//wrap
#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION
#include "bfm.h"
#include "Tools/Bagel/bfm_storage.hpp"


using namespace std;

int Test_Solver_BFM::run(){
  CCIO::cout << ".::: Test Solver BFM Moebius" 
	     <<std::endl;


  Mapping::init_shiftField();

  //Parameters
  int Nvol  =  CommonPrms::instance()->Nvol();
  double mq = 1.0;
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


  Dirac_optimalDomainWall_EvenOdd DWF_EO(b,c, -M5, mq, omega, &Dw_eo, &Dw_oe, &(conf_.data));
  
  /********************************************************

   * Setup DWF operator

   ********************************************************
   */

  bfmarg  dwfa;
  dwfa.node_latt[0]  = CommonPrms::instance()->Nx();
  dwfa.node_latt[1]  = CommonPrms::instance()->Ny();
  dwfa.node_latt[2]  = CommonPrms::instance()->Nz();
  dwfa.node_latt[3]  = CommonPrms::instance()->Nt();
  dwfa.verbose = 1;

  for(int mu=0;mu<4;mu++){
    if ( (CommonPrms::instance()->NPE(mu))>1 ) {
      dwfa.local_comm[mu] = 0;
    } else {
      dwfa.local_comm[mu] = 1;
    }
  }

  int threads =1;

  bfmarg::Threads(threads);
  bfmarg::UseCGdiagonalMee(1);

  dwfa.ScaledShamirCayleyTanh(mq,M5,Ls,ht_scale);
  dwfa.rb_precondition_cb=Even;
  dwfa.max_iter=10000;
  dwfa.residual=1.0e-8;
  bfm_dp linop;
  linop.init(dwfa);

  BFM_Storage BFM_interface(linop);
  BFM_interface.GaugeExport_to_BFM(conf_);
  
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
  // Tests the import/export routines
  ////////////////////////////////////////////////////////////////////////////////////
  Fermion_t psi_h[2];
  Fermion_t chi_h[2];
  Fermion_t tmp = linop.allocFermion();
  CCIO::cout << "tmp allocated \n";
  FermionField Exported(Nvol5d);
  for(int cb=0;cb<2;cb++){
    psi_h[cb] = linop.allocFermion();//half vector
    CCIO::cout << "phi_h["<<cb<<"] allocated \n";
    chi_h[cb] = linop.allocFermion();
    CCIO::cout << "chi_h["<<cb<<"] allocated \n";
    // Export fermion field to BFM
    for (int s = 0; s< Ls; s++){
      CCIO::cout << "Exporting Fermion to BFM cb = "<< cb <<" s = "<<s <<"\n";
      BFM_interface.FermionExport_to_BFM_5D(EO_source,psi_h[cb],cb,s);//third argument 0=even, 1=odd
    }
    // Import back to IroIro
    for (int s = 0; s< Ls; s++)
      BFM_interface.FermionImport_from_BFM_5D(Exported, psi_h[cb], cb,s);
  }
  //////////////////////////////////////////////////////////////////////////////////////
  // Check correctness of Import/Export routines
  FermionField Difference = EO_source;
  Difference -= Exported;

  CCIO::cout << "Difference source-exported = "<< Difference.norm() << "\n";
  CCIO::cout << "Norm Exp = "<< Exported.norm() << "\n";
  CCIO::cout << "Norm Src = "<< EO_source.norm() << "\n";
  ///////////////////////////////////////////////////////////////////////////////////////
  // Running the dslash
  int dag=0;
  int donrm=0;
  int cb=0;//even checkerboard to match conventions IroIro<->BFM
  //linop.MprecTilde(psi_h[Even],chi_h[Even],tmp,dag,donrm) ;

  //linop.Meo(psi_h[1-cb],chi_h[cb],cb,dag);
  //linop.Meo(psi_h[cb],chi_h[1-cb],1-cb,dag);

  //linop.Mooee(psi_h[cb],chi_h[cb],dag);
  //linop.Mooee(psi_h[1-cb],chi_h[1-cb],dag);

  //linop.MooeeInv(psi_h[cb],chi_h[cb],dag);
  //linop.MooeeInv(psi_h[1-cb],chi_h[1-cb],dag);
  
#pragma omp parallel
  {
#pragma omp for 
    for (int t=0;t<threads;t++){
      linop.CGNE_prec(chi_h[Even],psi_h[Even]);
    }
  }

  FermionField BFMsolution(Nvol5d);
  int vect4d_hsize = fe.size()/Ls;   
  for(int cb=0;cb<2;cb++){
    for (int s = 0; s< Ls; s++)
      BFM_interface.FermionImport_from_BFM_5D(BFMsolution, chi_h[cb], cb,s);
  }
  for (int s =0 ; s< Ls; s++){
    for (int i = 0; i < vect4d_hsize; i++){
      //      fe.data.set(i+2*s*vect4d_hsize,IroIroSol_even[i+vect4d_hsize*s]);
      fe.set(i+vect4d_hsize*s, BFMsolution.data[i+2*s*vect4d_hsize]);
    }
  }




  //Apply operator back on the inverse
  FermionField IroIroFull(Nvol5d);
  
  //Field IroIroSol_odd  = DWF_EO.mult_ee(DWF_EO.mult_oe(fe));
  //Field IroIroSol_even = DWF_EO.mult_oo(DWF_EO.mult_eo(fo));

  //Field IroIroSol_odd  = 0;
  Field IroIroSol_even = DWF_EO.mult_dag(DWF_EO.mult(fe));


  //////////////////////////////////////////////////////////////////////////////////  
  
  CCIO::cout << "IroIroSol_even.size()/Ls : "<< IroIroSol_even.size()/Ls << "\n";

  for (int s =0 ; s< Ls; s++){
    for (int i = 0; i < vect4d_hsize; i++){
      IroIroFull.data.set(i+2*s*vect4d_hsize,IroIroSol_even[i+vect4d_hsize*s]);
      //IroIroFull.data.set(i+(2*s+1)*vect4d_hsize, IroIroSol_odd[i+vect4d_hsize*s]);
      IroIroFull.data.set(i+(2*s+1)*vect4d_hsize, fo[i+vect4d_hsize*s]);
    }
  }
  


  for (int i = 0; i < IroIroFull.size() ; i++){
    double diff = abs(IroIroFull.data[i]-EO_source.data[i]);
    if (diff>1e-8) CCIO::cout << "*";
    CCIO::cout << "["<<i<<"] "<<IroIroFull.data[i] << "  "<<EO_source.data[i]
               << "  "<< diff << "\n";
  }
  
  Difference = EO_source;
  Difference -= IroIroFull;
  CCIO::cout << "Operator Difference BFM-IroIro = "<< Difference.norm() << "\n";
  
  return 0;
}
