/* @file test_PolyakovLoop.cpp
 * @brief implementation of the Test_PolyakovLoop.cpp class
 */
#include "Measurements/GaugeM/polyakovLoop.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include "test_PolyakovLoop.hpp"
#include <complex>
#include <cstring>
#include <iomanip>

int Test_PolyakovLoop::run(){
  XML::node pnode = input_.node;
  XML::descend(pnode,"PolyakovLoop");
  const char* dir_name = pnode.attribute("dir").value();

  Staples Staple;
  CCIO::cout<< "Plaquette : "<< Staple.plaquette(*(input_.gconf))<< std::endl;


  site_dir dir;
  if(     !strcmp(dir_name,"X")) dir = XDIR;
  else if(!strcmp(dir_name,"Y")) dir = YDIR;
  else if(!strcmp(dir_name,"Z")) dir = ZDIR;
  else if(!strcmp(dir_name,"T")) dir = TDIR;
  else {
    CCIO::cout<<"No valid direction available with name "<< dir_name << "\n";
    abort();
  }
  PolyakovLoop plp(dir);

  std::complex<double> pf = plp.calc_SUN(*(input_.gconf));
  double          pa = plp.calc_SUNadj(*(input_.gconf));
  
  CCIO::cout<<"Polyakov Loop ("<<dir << ") [fundamental representation]: ";
  CCIO::cout<< std::setw(20)<<pf.real()<<" "<< std::setw(20)<<pf.imag()<<"\n";
  CCIO::cout<<"Polyakov Loop [adjoint representation]:     ";
  CCIO::cout<< std:: setw(20) << pa <<"\n";
  
  return 0;
}
