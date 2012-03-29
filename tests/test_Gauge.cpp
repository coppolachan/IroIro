/*!
  @file test_Gauge.cpp
  @brief File for testing the Staples and Mapper classes
 */
#include "test_Gauge.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include "lib/Main/Geometry/mapping.hpp"

using namespace std;

int Test_Gauge::run() {
  CCIO::header("Test Gauge");

  Mapping::init_shiftField();

  int res1 = plaquette();
  int res2 = map_test();
  return (res1 && res2);
}

int Test_Gauge::plaquette(){
  CCIO::cout<<" Testing plaquette\n";
  Staples wl;
  CCIO::cout<<" Plaquette = "<<  wl.plaquette(d_conf) << endl;
  return 0;
}

int Test_Gauge::map_test(){
  using namespace Mapping;
  CCIO::cout<<" Testing maps\n";

  double loop_timer;
  Staples wl;
  int nid = Communicator::instance()->nodeid();

  GaugeField test_u;
  test_u = d_conf;               //create from existing data
  GaugeField1D test_u_1D(200);   //passing the lattice dimension
  FermionField test_v;

  CCIO::cout<<"\n Testing common field templates\n Local sizes:\n";
  CCIO::cout<<" GaugeFieldType size   : "<< test_u.data.size()   <<"\n";
  CCIO::cout<<" GaugeField1DType size : "<< test_u_1D.data.size()<<"\n";
  CCIO::cout<<" FermionField size     : "<< test_v.data.size()   <<"\n";

  CCIO::cout<<" plaq (original) = "
	    <<wl.plaquette(test_u)<< endl; 
  CCIO::cout<<" plaq (+x map  ) = "
	    <<wl.plaquette(shiftField(test_u,0,Forward()))<< endl; 
  CCIO::cout<<" plaq (-x map  ) = "
	    <<wl.plaquette(shiftField(test_u,0,Backward()))<< endl; 
  CCIO::cout<<" plaq (+y map  ) = "
	    <<wl.plaquette(shiftField(test_u,1,Forward()))<< endl; 
  CCIO::cout<<" plaq (-y map  ) = "
	    <<wl.plaquette(shiftField(test_u,1,Backward()))<< endl; 
  CCIO::cout<<" plaq (+z map  ) = "
	    <<wl.plaquette(shiftField(test_u,2,Forward()))<< endl; 
  CCIO::cout<<" plaq (-z map  ) = "
	    <<wl.plaquette(shiftField(test_u,2,Backward()))<< endl; 
  CCIO::cout<<" plaq (+t map  ) = "
	    <<wl.plaquette(shiftField(test_u,3,Forward()))<< endl; 
  CCIO::cout<<" plaq (-t map  ) = "
	    <<wl.plaquette(shiftField(test_u,3,Backward()))<< endl; 
  

  return 0;
}

