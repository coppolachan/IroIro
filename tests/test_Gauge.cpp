/*!
  @file test_Gauge.cpp
   
  @brief File for testing the gauge Field

 */

#include "test_Gauge.hpp"
#include "Measurements/GaugeM/staples.h"
#include "Communicator/comm_io.hpp"
#include "include/common_fields.hpp"


using namespace std;
using namespace Format;

int Test_Gauge::run() {
  CCIO::header("Test Gauge");

  int res1 = plaquette();
  int res2 = shift();
  return (res1 && res2);

}

int Test_Gauge::plaquette(){
  Staples wl(d_conf.get_Format());
  double plaq = wl.plaquette(d_conf.data);

  CCIO::cout<<" Plaq = "<<plaq<<endl;
  return 0;
}

int Test_Gauge::shift(){
  CCIO::cout<<" Testing shifter:\n";
  Format_G FormatObj = d_conf.get_Format();
  ShiftField_up<Format_G> upt(d_conf.data,&FormatObj,3); //+t
  ShiftField_dn<Format_G> umx(d_conf.data,&FormatObj,0); //-x

  Staples wl(FormatObj);
  int nid = Communicator::instance()->nodeid();

  CCIO::cout <<" plaq (original) = "<<wl.plaquette(d_conf.data)<<" nodeid= "<<nid << endl;
  CCIO::cout <<" plaq (+t sfted) = "<<wl.plaquette(upt)     <<" nodeid= "<<nid << endl;
  CCIO::cout <<" plaq (-x sfted) = "<<wl.plaquette(umx)     <<" nodeid= "<<nid << endl;


  GaugeFieldType test_u;
  test_u = d_conf;   //create from existing data
  GaugeField1DType test_u_1D(200);   //passing the lattice dimension
  FermionField test_v;
  //FermionFieldExtraDim test_v_5D; // should not compile
  FermionFieldExtraDim test_v_5D(8); 


  CCIO::cout <<"\n Testing common field templates\n";
  CCIO::cout <<" GaugeFieldType size   : " << test_u.data.size() << "\n";
  CCIO::cout <<" GaugeField1DType size : " << test_u_1D.data.size() << "\n";
  CCIO::cout <<" FermionField size     : " << test_v.data.size() << "\n";
  CCIO::cout <<" FermionFieldExtraDim  : " << test_v_5D.data.size() << "\n";


  return 0;
}

