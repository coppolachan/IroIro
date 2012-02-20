/*!
  @file test_Gauge.cpp
   
  @brief File for testing the gauge Field

 */

#include "test_Gauge.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include "Communicator/comm_io.hpp"
#include "include/common_fields.hpp"

#include "Main/Geometry/mapper.hpp"


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
  double plaq = wl.plaquette(d_conf);

  CCIO::cout<<" Plaquette = "<<plaq<<endl;
  return 0;
}

int Test_Gauge::shift(){
  CCIO::cout<<" Testing shifter:\n";
  Format_G FormatObj = d_conf.get_Format();
  ShiftField_up<Format_G> upt(d_conf.data,&FormatObj,3); //+t
  ShiftField_dn<Format_G> umx(d_conf.data,&FormatObj,0); //-x
  Mapper shift;

  double loop_timer;
  TIMING_START;
  for (int loop = 0; loop< 1000; ++loop) {
    ShiftField_up<Format_G> stest(d_conf.data,&FormatObj,2); //+z
    std::valarray<double> va_stest = upt.getva();
  }
  TIMING_END(loop_timer);
  CCIO::cout <<"ShiftField Loop timing: "<< loop_timer<< endl;

  Staples wl(FormatObj);
  int nid = Communicator::instance()->nodeid();
 
  CCIO::cout <<" plaq (original) = "<<wl.plaquette(d_conf.data)<<" nodeid= "<<nid << endl;
  CCIO::cout <<" plaq (+t sfted) = "<<wl.plaquette(upt)     <<" nodeid= "<<nid << endl;
  CCIO::cout <<" plaq (-x sfted) = "<<wl.plaquette(umx)     <<" nodeid= "<<nid << endl;

  gettimeofday(&start,NULL);
  for (int loop = 0; loop< 100; ++loop) 
    wl.upper(d_conf.data,0,1);
  gettimeofday(&end,NULL);
  loop_timer = (end.tv_sec - start.tv_sec) * 1000.0;                          
  loop_timer = loop_timer + (end.tv_usec - start.tv_usec) / 1000.0 ;  // us to ms
    
  CCIO::cout <<" Old staple Loop timing: "<< loop_timer<< endl;

  CCIO::cout <<" lower norm     = "<<(wl.upper(d_conf.data,0,1)).norm() <<" nodeid= "<<nid << endl;

  GaugeFieldType test_u;
  test_u = d_conf;                   //create from existing data
  GaugeField1DType test_u_1D(200);   //passing the lattice dimension
  FermionField test_v;
  //FermionFieldExtraDim test_v_5D; // should not compile
  FermionFieldExtraDim test_v_5D(8); 


  CCIO::cout <<"\n Testing common field templates\n";
  CCIO::cout <<" GaugeFieldType size   : " << test_u.data.size() << "\n";
  CCIO::cout <<" GaugeField1DType size : " << test_u_1D.data.size() << "\n";
  CCIO::cout <<" FermionField size     : " << test_v.data.size() << "\n";
  CCIO::cout <<" FermionFieldExtraDim  : " << test_v_5D.data.size() << "\n";

  gettimeofday(&start,NULL);
  for (int loop = 0; loop< 1000; ++loop) 
    GaugeFieldType map_test = shift(test_u, 0,  Forward);
  gettimeofday(&end,NULL);
  loop_timer = (end.tv_sec - start.tv_sec) * 1000.0;                          
  loop_timer = loop_timer + (end.tv_usec - start.tv_usec) / 1000.0 ;  // us to ms
    
  CCIO::cout <<"Mapper Loop timing: "<< loop_timer<< endl;

  
  CCIO::cout <<" plaq (original) = "<<wl.plaquette(test_u)<< endl; 
  CCIO::cout <<" plaq (+x map  ) = "<<wl.plaquette(shift(test_u, 0,  Forward))<< endl; 
  
  CCIO::cout <<" plaq (-x map  ) = "<<wl.plaquette(shift(test_u, 0, Backward))<< endl; 
  CCIO::cout <<" plaq (+y map  ) = "<<wl.plaquette(shift(test_u, 1,  Forward))<< endl; 
  CCIO::cout <<" plaq (-y map  ) = "<<wl.plaquette(shift(test_u, 1, Backward))<< endl; 
  CCIO::cout <<" plaq (+z map  ) = "<<wl.plaquette(shift(test_u, 2,  Forward))<< endl; 
  CCIO::cout <<" plaq (-z map  ) = "<<wl.plaquette(shift(test_u, 2, Backward))<< endl; 
  CCIO::cout <<" plaq (+t map  ) = "<<wl.plaquette(shift(test_u, 3,  Forward))<< endl; 
  CCIO::cout <<" plaq (-t map  ) = "<<wl.plaquette(shift(test_u, 3, Backward))<< endl; 
  

  // Specialized Staples
  gettimeofday(&start,NULL);
  for (int loop = 0; loop< 100; ++loop) 
    wl.upper(d_conf,0,1);
  gettimeofday(&end,NULL);
  loop_timer = (end.tv_sec - start.tv_sec) * 1000.0;                          
  loop_timer = loop_timer + (end.tv_usec - start.tv_usec) / 1000.0 ;  // us to ms
  CCIO::cout <<"New staple lower Loop timing: "<< loop_timer<< endl;

  CCIO::cout <<" Map staple norm     = "<<(wl.upper(d_conf,0,1)).norm() <<" nodeid= "<<nid << endl;
  
  for(int i=0;i<NDIM_-1;++i){
    int j = (i+1)%(NDIM_-1);  
    Field diff_test = wl.upper(d_conf,i,j).data;
    diff_test -= wl.upper(d_conf.data,i,j);
    CCIO::cout <<" Test - Map diff norm    = "<< diff_test.norm() <<endl;
  }
  return 0;
}

