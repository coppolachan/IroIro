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
  /*
  /////////////////////////////
  XML::node QuarkProp_node = DWFnode;
  XML::descend(QuarkProp_node, "QuarkPropagator");
  QPropDWFFactory  QP_DomainWallFact(QuarkProp_node);//uses specific factory (this is a test program specific for DWF)
  QpropDWF* QuarkPropDW = static_cast<QpropDWF*>(QP_DomainWallFact.getQuarkProp(conf_));
  // the prevoius static_cast is absolutely safe since we know exaclty what class we are creating
  
  ///////////////////////////////////////
  */
  CCIO::cout << ".::: Test Solver BFM" 
	     <<std::endl;


  Mapping::init_shiftField();
  int Nvol =  CommonPrms::instance()->Nvol();
  // Random number generator initialization
  unsigned long init[4]={0x123, 0x234, 0x345, 0x456};
  int length=4;
  RandNum_MT19937 rand(init, length);
  
  vector<int> spos(4,0);
  spos[0] = 0;
  spos[1] = 0;
  spos[2] = 0;
  spos[3] = 0;

  double alpha = 5.0;
  Source_local<Format::Format_F> src(spos,CommonPrms::instance()->Nvol());
  //Source_Gauss<Format::Format_F> src(spos,alpha, Nvol);
  
  prop_t sq;

  //  QuarkPropDW->calc(sq,src);
  
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
  dwfa.time_report_iter=100;


  for(int mu=0;mu<4;mu++){
    if ( (CommonPrms::instance()->NPE(mu))>1 ) {
      dwfa.local_comm[mu] = 0;
    } else {
      dwfa.local_comm[mu] = 1;
    }
  }
  double mq = 0.0;
  double M5 = 1.6;
  int Ls = 8;
  double ht_scale = 2.0;

  //  dwfa.ScaledShamirCayleyTanh(mq,M5,Ls,ht_scale);
  dwfa.pWilson(mq);

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
  Field fe = src.mksrc(esec,2,0);
  Field fo = src.mksrc(osec,2,0);
#ifdef DEBUG
  valarray<double> phi(source.format.size()/2);
  MPrand::mp_get_gauss(phi, rand);
  Field fe(phi);
  Field fo(phi);
#endif
  FermionField EO_source;
  
  CCIO::cout << "fe size: "<< fe.size() <<"\n";
  //copy fields in the eo source CB blocked
  assert(fe.size() == fo.size());
  for (int i = 0; i < fe.size(); i++){
    EO_source.data.set(i, fe[i]);
    //EO_source.data.set(i+fe.size(), fo[i]);
  }

  // Tests the import/export routines
  ///////////////////////////////////////
  Fermion_t psi_h;
  Fermion_t chi_h;
  psi_h = linop.allocFermion();//half vector
  chi_h = linop.allocFermion();
  // Export fermion field to BFM
  BFM_interface.FermionExport_to_BFM(EO_source,psi_h,0);//third argument is the CB 0=even, 1=odd
 
  // Import back to IroIro
  FermionField Exported;
  BFM_interface.FermionImport_from_BFM(Exported, psi_h, 0);

  // Check correctness of Import/Export routines
  FermionField Difference = EO_source;
  Difference -= Exported;

  CCIO::cout << "Difference source-exported = "<< Difference.norm() << "\n";

  // Running the dslash
  for(int dag=0;dag<1;dag++){
    linop.dslash(psi_h,
		 chi_h,
		 1, dag);
  }
   
  FermionField BFMsolution; 
  BFM_interface.FermionImport_from_BFM(BFMsolution, chi_h, 0);

  Dirac_Wilson_EvenOdd Deo(mq, &(conf_.data));
  Dirac_Wilson Kernel(mq, &(conf_.data));
  Field IroIroSol_eo = Deo.mult_eo(fo);
  Field IroIroSol_oe = Deo.mult_oe(fe);

  //Normalization conventions
  IroIroSol_eo *= (-1.0/Deo.getKappa());
  IroIroSol_oe *= (-1.0/Deo.getKappa());

  FermionField IroIroSol;

  for (int i = 0; i < IroIroSol_oe.size(); i++){
    IroIroSol.data.set(i, IroIroSol_oe[i]);
    //EO_source.data.set(i+fe.size(), fo[i]);
#ifdef DEBUG
    double diff = abs(IroIroSol.data[i]-BFMsolution.data[i]);
    if (diff>1e-8) CCIO::cout << "*";
    CCIO::cout << "["<<i<<"] "<<IroIroSol.data[i] << "  "<<BFMsolution.data[i] 
	       << "  "<< diff << "\n";
#endif
  }

  Difference = BFMsolution;
  Difference -= IroIroSol;
  CCIO::cout << "Operator Difference BFM-IroIro_oe = "<< Difference.norm() << "\n";
  
  return 0;
}
