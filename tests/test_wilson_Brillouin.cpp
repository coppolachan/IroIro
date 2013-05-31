/*!
 * @file test_wilson_Brillouin.cpp
 * @brief Tests for the Dirac_Wilson_Brillouin class
 */
#include "test_wilson_Brillouin.hpp"
#include "Measurements/FermionicM/fermion_meas_factory_abs.hpp"

#include "Dirac_ops/dirac_wilson_Brillouin.hpp"
#include <stdio.h>
#include <time.h>

using namespace std;
using namespace Format;

int Test_Wilson_Brillouin::run(){

  Dirac_Wilson_Brillouin D(0.0,&(conf_.data));
  
  Field f(D.fsize());
  int Nvol = CommonPrms::instance()->Nvol();
  ffmt_t ff(Nvol);
  
  f.set(0,1.0);
  
  Field w = D.mult(f);
  
  for(int site=0;site<Nvol;++site){
    int x  = SiteIndex::instance()->c_x(site);
    int y  = SiteIndex::instance()->c_y(site);
    int z  = SiteIndex::instance()->c_z(site);
    int t  = SiteIndex::instance()->c_t(site);
   
    cout<<"("<<x<<","<<y<<","<<z<<","<<t<<")="<<w[ff.index(0,site)]<<endl;
  }
}
