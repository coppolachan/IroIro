/*!
 * @file test_Solver_BFM.cpp
 * @brief Definition of classes for testing the BFM classes
 */
#include "test_Solver_BFM.hpp"
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


  vector<int> spos(4,0); 
  Source_local<Format::Format_F> src(spos,CommonPrms::instance()->Nvol());
  
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


  for(int mu=0;mu<4;mu++){
    if ( (CommonPrms::instance()->NPE(mu))>1 ) {
      dwfa.local_comm[mu] = 0;
    } else {
      dwfa.local_comm[mu] = 1;
    }
  }
  double mq = 0.1;
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

  // Tests the import/export routines
  ///////////////////////////////////////
  Fermion_t psi_h;
  Fermion_t chi_h;
  psi_h = linop.allocFermion();
  chi_h = linop.allocFermion();
  // Export fermion field to BFM
  BFM_interface.FermionExport_to_BFM(source,psi_h,0);
 
  // Import back to IroIro
  FermionField Exported;
  BFM_interface.FermionImport_from_BFM(Exported, psi_h, 0);

  // Check correctness of Import/Export routines
  FermionField Difference = source;
  Difference -= Exported;

  CCIO::cout << "Difference source-exported = "<< Difference.norm() << "\n";

  // Running the dslash
  for(int dag=0;dag<2;dag++){
    linop.dslash(psi_h,
		 chi_h,
		 1, dag);
  }

  return 0;
}
