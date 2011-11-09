//--------------------------------------------------------------------
// mdExec_leapfrog.cpp
//--------------------------------------------------------------------
#include "mdExec_leapfrog.h"
#include "include/field.h"
#include "Tools/sunMat.h"
#include "Measurements/GaugeM/staples.h"
#include <iomanip>
using namespace std;
using namespace Format;

inline SUNmat MDexec_leapfrog::
u(const Field& g,int site,int dir)const{
  return SUNmat(g[gf_.cslice(0,site,dir)]);
}

void MDexec_leapfrog::
update_U(const Field& P,double ep)const{
  
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
      // au *= u(U,site,m);
      // U.set(gf_.islice(site,m),au.reunit().getva());
      au *= u(*CommonField,site,m);
      CommonField->set(gf_.islice(site,m),au.reunit().getva());

    }
  }
  //if(stpl_) printf("plaq1=%.12f\n",stpl_->plaquette(U));    
}

void MDexec_leapfrog::
update_P(Field& P,int lv,double ep) const{
  for(int a = 0; a < as_[lv].size(); ++a){
    Field fse = as_[lv].at(a)->md_force();
    //cout<<"fse.norm()="<<fse.norm()<<endl;
    fse *= ep;
    P -= fse; 
  }
}

void MDexec_leapfrog::
init(vector<int>& clock,Field& P,const Field& U,
     const RandNum& rand) const{

  //Initialize the Common field to U
  *CommonField = U;

  int nodeid=Communicator::instance()->nodeid();
  clock.resize(as_.size(),0.0);

  for(int lv = 0; lv< as_.size(); ++lv){
    for(int id = 0; id < as_.at(lv).size(); ++id){
      
      if(nodeid==0) 
	cout<<"initializing MD steps level= "<< lv <<" id= "<< id<<endl;
      as_[lv].at(id)->init(P,rand);
    }
  }
}

double MDexec_leapfrog::calc_H()const{

  int nodeid=Communicator::instance()->nodeid();
  double H=0.0;
  for(int lv = 0; lv< as_.size(); ++lv)
    for(int id = 0; id < as_.at(lv).size(); ++id){
      if(nodeid==0) cout<<" id="<<id<<" lv="<<lv<<endl;
      H+= as_[lv].at(id)->calc_H();
    }
  return H;
}




void MDexec_leapfrog::
integrator_step(Field& P,
		int cl,
		std::vector<int>& clock,
		int nodeid) const{
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
      update_P(P,cl,eps/2);
      ++clock[cl];
      if(nodeid==0){
	for(int l=0; l<=cl;++l) cout<<"   ";
	cout<<"P "<<static_cast<double>(clock[cl])/2 <<endl;
      }
    }
    if(cl == fl){          // lowest level 
      update_U(P,eps);
      if(nodeid==0){
	for(int l=0; l<=cl;++l) cout<<"   ";
	cout<<"U "<< static_cast<double>(clock[cl]+1)/2 <<endl;
      }
    }else{                 // recursive function call 
      integrator_step(P,cl+1,clock,nodeid);
    }
    if(clock[cl] == fin){  // final half step
      update_P(P,cl,eps/2);
      
      ++clock[cl];
      if(nodeid==0){
	for(int l=0; l<=cl;++l) cout<<"   ";
	cout<<"P "<< static_cast<double>(clock[cl])/2 <<endl;
      }
    }else{                  // bulk step
      update_P(P,cl,eps);
      
      clock[cl]+=2;
      if(nodeid==0){
	for(int l=0; l<=cl;++l) cout<<"   ";
	cout<<"P "<< static_cast<double>(clock[cl])/2 <<endl;
      }
    }
  }
  
}


void MDexec_leapfrog::
integrator(Field& P,int cl,std::vector<int>& clock) const{
  // cl  : current level
  // fl  : final level
  // eps : current step size
  
  int nodeid=Communicator::instance()->nodeid();
  
  for(int step=0; step< Params.MDsteps; ++step){   // MD step 
    if(nodeid==0) cout<<"MDstep = "<< step << endl;
    integrator_step(P,cl,clock,nodeid);
  }
  
};
