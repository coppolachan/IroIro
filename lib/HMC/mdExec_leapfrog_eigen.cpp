//--------------------------------------------------------------------
// mdExec_leapfrog_eigen.cpp
//--------------------------------------------------------------------
#include "mdExec_leapfrog_eigen.h"
#include "include/field.h"
#include "Tools/sunMat.h"
#include "Measurements/GaugeM/staples.h"
#include "EigenModes/eigenModes_IRL.h"
#include <iomanip>

using namespace std;
using namespace Format;

inline SUNmat MDexec_leapfrog_eigen::u(const Field& g,
				       int site,
				       int dir) const{
  return SUNmat(g[gf_.cslice(0,site,dir)]);
}

void MDexec_leapfrog_eigen::update_U(Field& P,
				     double ep) const{
  
  using namespace SUNmat_utils;
  
  //if(stpl_) printf("plaq0=%.12f\n",stpl_->plaquette(U));    
  
  int Ndim = CommonPrms::instance()->Ndim();
  int Nvol = CommonPrms::instance()->Nvol();
  
  const SUNmat I = unity();
  
  for(int m = 0; m < Ndim; ++m){
    for(int site=0; site<Nvol; ++site){
      SUNmat au = I;
      for(int k = Params.Nexp; k > 0; --k){
        //cout << k << ", " << ep << ", " << ep/k << endl;
	au *= ep/k;
	au *= u(P,site,m);
	au += I;
      }
      au *= u(*CommonField,site,m);
      CommonField->set(gf_.islice(site,m),au.reunit().getva());
    }
  }
  //if(stpl_) printf("plaq1=%.12f\n",stpl_->plaquette(U));    
}

void MDexec_leapfrog_eigen::update_P(Field& P,
				     int lv,
				     double ep) const{

  EigenData ed;

  for(int a = 0; a < as_[lv].size(); ++a){
    if(es_[lv].at(a)) es_[lv].at(a)->calc(ed);

    Field fse = as_[lv].at(a)->md_force(&ed);
    //cout<<"fse.norm()="<<fse.norm()<<endl;
    fse *= ep;
    P -= fse; 
  }
}

void MDexec_leapfrog_eigen::
init(vector<int>& clock,Field& P,const Field& U,
     const RandNum& rand) const{

  *CommonField = U;

  clock.resize(as_.size(),0.0);
  EigenData ed;
  
  for(int lv = 0; lv< as_.size(); ++lv){
    for(int a = 0; a < as_[lv].size(); ++a){
      if(es_[lv].at(a)) es_[lv].at(a)->calc(ed,U);

      cout<<"initializing MD steps level= "<< lv <<" id= "<< a<<endl;
      as_[lv].at(a)->init(P,rand,U,&ed);
    }
  }
}

double MDexec_leapfrog_eigen::calc_H()const{
  double H;
  for(int lv = 0; lv< as_.size(); ++lv)
    for(int id = 0; id < as_.at(lv).size(); ++id)
      H+= as_[lv].at(id)->calc_H();
  return H;
}

void MDexec_leapfrog_eigen::
integrator(Field& P,Field& U,int cl,std::vector<int>& clock) const{
  // cl  : current level
  // fl  : final level
  // eps : current step size
  
  int fl = as_.size() -1;
  double eps = Params.step_size;
  for(int l=0; l<=cl; ++l) eps/= Nrel_[l];

  int fin = 1;
  for(int l=0; l<=cl; ++l) fin*= Nrel_[l];
  fin = 2*Params.MDsteps*fin -1;

  for(int e=0; e<Nrel_[cl]; ++e){

    if(clock[cl] == 0){    // initial half step 
      update_P(P,U,cl,eps/2);
      ++clock[cl];
      for(int l=0; l<=cl;++l) cout<<"   ";
      cout<<"P "<<static_cast<double>(clock[cl])/2.0 <<endl;
    }
    if(cl == fl){          // lowest level 
      update_U(U,P,eps);
      for(int l=0; l<=cl;++l) cout<<"   ";
      cout<<"U "<< static_cast<double>(clock[cl]+1)/2.0 <<endl;

    }else{                 // recursive function call 
      integrator(P,U,cl+1,clock);
    }
    if(clock[cl] == fin){  // final half step
      update_P(P,U,cl,eps/2);

      ++clock[cl];
      for(int l=0; l<=cl;++l) cout<<"   ";
      cout<<"P "<< static_cast<double>(clock[cl])/2.0 <<endl;

    }else{                  // bulk step
      update_P(P,U,cl,eps);

      clock[cl]+=2;
      for(int l=0; l<=cl;++l) cout<<"   ";
      cout<<"P "<< static_cast<double>(clock[cl])/2.0 <<endl;
    }
  }
}

