/*!
 * @file test_wilson_Brillouin.cpp
 * @brief Tests for the Dirac_Wilson_Brillouin class
 */
#include "test_wilson_Brillouin.hpp"
#include "Measurements/FermionicM/fermion_meas_factory_abs.hpp"

#include "Dirac_ops/dirac_wilson_Brillouin.hpp"
//#include "Solver/solver_CG.hpp"
//#include "Solver/solver_BiCGStab.hpp"
//#include "Measurements/GaugeM/staples.hpp"
//#include "Measurements/FermionicM/qprop.hpp"
//#include "Measurements/FermionicM/source_types.hpp"
//#include "Measurements/FermionicM/meson_correlator.hpp"
//#include "Tools/randNum_MT19937.h"
#include <stdio.h>
#include <time.h>

using namespace std;
using namespace Format;

/*
timespec diff(timespec start, timespec end){
  timespec temp;
  if ((end.tv_nsec-start.tv_nsec)<0) {
    temp.tv_sec = end.tv_sec-start.tv_sec-1;
    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec-start.tv_sec;
    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return temp;
}
*/

int Test_Wilson_Brillouin::run(){

  Dirac_Wilson_Brillouin D(0.0,&(conf_.data));
  
  Field f(D.fsize());
  Format::Format_F ff = D.get_fermionFormat();
  
  f.set(0,1.0);
  
  Field w = D.mult(f);
  
  int Nvol = CommonPrms::instance()->Nvol();
  
  for(int site=0;site<Nvol;++site){
    int x  = SiteIndex::instance()->c_x(site);
    int y  = SiteIndex::instance()->c_y(site);
    int z  = SiteIndex::instance()->c_z(site);
    int t  = SiteIndex::instance()->c_t(site);
   
    cout<<"("<<x<<","<<y<<","<<z<<","<<t<<")="<<w[ff.index(0,site)]<<endl;
  }
}
