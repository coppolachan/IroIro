/*!
 * @file test_wilson_Brillouin.cpp
 * @brief Tests for the Dirac_Wilson_Brillouin class
 Time-stamp: <2013-06-04 14:05:36 noaki>
 */
#include "test_wilson_Brillouin_Imp.hpp"
#include "Measurements/FermionicM/quark_prop_meas_factory.hpp"
#include "Dirac_ops/dirac_wilson_Brillouin.hpp"
#include "include/numerical_const.hpp"
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <fstream>
#include <cstdio>

using namespace std;
using namespace Format;

int Test_Wilson_Brillouin_Imp::run(){
  int Nvol = CommonPrms::instance()->Nvol();
  Dirac_Wilson_Brillouin D(0.0,&(conf_.data),Improved);
  Field f(D.fsize());
  fstream file; 
  Format::Format_F ff(Nvol);
  const int L = CommonPrms::instance()->Lx();
  double sum = 0.0;
  
  // set source vector
  f.set(0,1.0);

  // multiple Dirac operator
  Field w = D.mult(f);

  //check
  file.open("coefficient.d",ios::out);
  for(int site=0;site<Nvol;++site){
    int x  = SiteIndex::instance()->c_x(site);
    int y  = SiteIndex::instance()->c_y(site);
    int z  = SiteIndex::instance()->c_z(site);
    int t  = SiteIndex::instance()->c_t(site);
    
    int gx  = SiteIndex::instance()->global_x(x);
    int gy  = SiteIndex::instance()->global_y(y);
    int gz  = SiteIndex::instance()->global_z(z);
    int gt  = SiteIndex::instance()->global_t(t);
    
    if(w[ff.index(0,site)]!=0){
      file<<"("<<gx<<","<<gy<<","<<gz<<","<<gt<<")="<<w[ff.index(0,site)]<<endl;
    }
  }
   file.close();  

  /*-------------------------------Fourier Tr --------------------------------*/
  // set parameters
  const int P_Vol = 30;
  const double Stp = 0.1;
  const int V = P_Vol*P_Vol*P_Vol*P_Vol;
  double px,py,pz,pt;

  double Fr[V];
  double Fi[V];
  double R[V];

  //write parameters
  cout<<"Volume = "<<Nvol<<endl;
  cout<<"P_Vol = "<<P_Vol<<endl;
  cout<<"L = "<<L<<endl;

  // main part of DFT
  for(int it = 0; it<P_Vol; ++it){
    for(int iz = 0; iz<P_Vol; ++iz){
      for(int iy = 0; iy<P_Vol; ++iy){	    
	for(int ix = 0; ix<P_Vol; ++ix){
	  //	  cout<<"(px,py,pz,pt)=("<<ix<<","<<iy<<","<<iz<<","<<it<<")"<<endl;	  
	  // determine each directions momentum from the step size
	  px = ix*Stp*2.0*PI/L;
	  py = iy*Stp*2.0*PI/L;
	  //	  pz = iz*Stp*2.0*PI/L;
	  //pt = it*Stp*2.0*
	  pz = pt = 0.0;

	  for(int site=0;site<Nvol;++site){
	    int x  = SiteIndex::instance()->c_x(site);
	    int y  = SiteIndex::instance()->c_y(site);
	    int z  = SiteIndex::instance()->c_z(site);
	    int t  = SiteIndex::instance()->c_t(site);
	    
	    int gx  = SiteIndex::instance()->global_x(x);
	    int gy  = SiteIndex::instance()->global_y(y);
	    int gz  = SiteIndex::instance()->global_z(z);
	    int gt  = SiteIndex::instance()->global_t(t);
	    
	  // periodicity
	    if(gx < L/2){
	      gx = gx + L/2;
	    }else{
	      gx = gx - L/2 ;
	    }
	    if(gy < L/2){
	      gy = gy + L/2;
	    }else{
	      gy = gy - L/2;
	    }
	    if(gz < L/2){
	      gz = gz + L/2;
	    }else{
	      gz = gz - L/2;
	    }
	    if(gt < L/2){
	      gt = gt + L/2;
	    }else{
	      gt = gt - L/2;
	    }

	    double srd = px*gx+py*gy+pz*gz+pt*gt;
	    unsigned long int p = ix + P_Vol*iy + P_Vol*P_Vol*iz + P_Vol*P_Vol*P_Vol*it;	    

	    Fr[p] += cos(srd)*w[ff.index(0,site)];
	    Fi[p] += -sin(srd)*w[ff.index(0,site)];
	    
	  }
	}
      }
    }
  }

  // Brillouin operator in momentum space
  file.open("function.d",ios::out);
  for(int it = 0; it<P_Vol; ++it){
    for(int iz = 0; iz<P_Vol; ++iz){
      for(int iy = 0; iy<P_Vol; ++iy){	    
	for(int ix = 0; ix<P_Vol; ++ix){
	  px = ix*Stp*2.0*PI/L;
	  py = iy*Stp*2.0*PI/L;
	  //pz = iz*Stp*2.0*PI/L;
	  //pt = it*Stp*2.0*PI/L;
	  pz = pt = 0.0;

	  unsigned long int p = ix + P_Vol*iy + P_Vol*P_Vol*iz + P_Vol*P_Vol*P_Vol*it;	    

	  R[p] = sin(px)*(cos(py)+2.0)*(cos(pz)+2.0)*(cos(pt)+2.0)/27 ; // x-derivative 
	  //R[p] = 4.0*pow(cos(px/2)*cos(py/2)*cos(pz/2)*cos(pt/2),2.0)-4.0; //laplacian
	  file<<px<<"  "<<py<<"  "<<pz<<"  "<<pt<<"  "<<R[p]<<endl;
	}
      }
    }
  }
  file.close();
  
  // write datas
  file.open("FourierForm.d",ios::out);
  for(int it = 0;it<P_Vol;++it){
    for(int iz = 0;iz<P_Vol;++iz){
      for(int iy = 0;iy<P_Vol;++iy){	    
	for(int ix = 0;ix<P_Vol;++ix){

	  unsigned long int p = ix + P_Vol*iy + P_Vol*P_Vol*iz + P_Vol*P_Vol*P_Vol*it;	    

	  double FF = sqrt(pow(Fr[p],2.0)+pow(Fi[p],2.0));
	  double RR = sqrt(pow(R[p],2.0));  
	  printf("p = %lu RR= %f  FF= %f\n",p,RR,FF);
	  file<<ix<<"  "<<iy<<"  "<<RR<<"  "<<FF<<endl;
	}
      }
    }
  }
  file.close();
}
