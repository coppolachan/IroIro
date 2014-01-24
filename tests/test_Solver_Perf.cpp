/*!
 * @file test_Solver_Perf.cpp
 * @brief Definition of classes for testing the Dirac_DomainWall classes and factories
 */
#include "test_Solver_Perf.hpp"
#include "Measurements/FermionicM/quark_prop_meas_factory.hpp"
#include "Measurements/FermionicM/qprop_DomainWall.hpp"
#include "Measurements/FermionicM/source_types.hpp"

using namespace std;

int Test_Solver_Performance::run(){
  InputConfig config(&conf_);
  /////////////////////////////
  XML::node QuarkProp_node = DWFnode;
  XML::descend(QuarkProp_node, "QuarkPropagator");
  QPropDWFFactory  QP_DomainWallFact(DWFnode);//uses specific factory (this is a test program specific for DWF)
  QpropDWF* QuarkPropDW 
    = static_cast<QpropDWF*>(QP_DomainWallFact.getQuarkProp(config));
  // the prevoius static_cast is absolutely safe since we know exaclty what class we are creating
  
  ///////////////////////////////////////
  CCIO::cout << ".::: Test Solver Performance" 
	     <<std::endl;

  vector<int> spos(4,0); 
  Source_local<Format::Format_F> src(spos,CommonPrms::instance()->Nvol());
  
  prop_t sq;
  QuarkPropDW->calc(sq,src);
  
  return 0;
}
