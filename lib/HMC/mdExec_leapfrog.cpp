//--------------------------------------------------------------------
// mdExec_leapfrog.cpp
//--------------------------------------------------------------------
#include "mdExec_leapfrog.h"
#include "include/field.h"
#include "Tools/sunMat.h"
#include <iomanip>

using namespace std;
using namespace Format;

inline SUNmat MDexec_leapfrog::
u(const Field& g,int site,int dir)const{
  return SUNmat(g[gf_.cslice(0,site,dir)]);
}

void MDexec_leapfrog::update_U(double ep){
  using namespace SUNmat_utils;

  int Ndim = CommonPrms::instance()->Ndim();
  int Nvol = CommonPrms::instance()->Nvol();

  const SUNmat I = unity();
  
  for(int m = 0; m < Ndim; ++m){
    for(int site=0; site<Nvol; ++site){
      SUNmat au = I;
      for(int k = Params.Nexp; k > 0; --k){
        //cout << k << ", " << ep << ", " << ep/k << endl;
	au *= ep/k;
	au *= u(P_,site,m);
	au += I;
      }
      au *= u(*U_,site,m);
      U_->set(gf_.islice(site,m),au.reunit().getva());
    }
  }
}

void MDexec_leapfrog::update_P(int lv,double ep){
  for(int a = 0; a < as_[lv].size(); ++a){
    Field fce = as_[lv].at(a)->md_force();
    fce *= ep;
    P_-= fce; 
  }
}

void MDexec_leapfrog::
init(vector<int>& clock,const Field& U,const RandNum& rand){
  clock.resize(as_.size(),0.0);  

  *U_= U;                       // initialize U_ (common to actions) to U
  MDutils::md_mom(P_,rand,gf_); // initialize P_ 

  double pnorm = P_.norm();
  CCIO::cout<<"P_->norm()= "<<pnorm<<endl;

  for(int lv = 0; lv< as_.size(); ++lv){
    for(int id = 0; id < as_.at(lv).size(); ++id){
      CCIO::cout<<"initializing MD steps level= "<< lv <<" id= "<< id<<endl;
      as_[lv].at(id)->init(rand);
    }
  }
}

double MDexec_leapfrog::calc_H()const{
  using namespace SUNmat_utils;
  // kinetic term
  double H_local = 0.0;
  Format::Format_G gf(CommonPrms::instance()->Nvol());

  for(int site = 0; site < CommonPrms::instance()->Nvol(); ++site){
    for(int dir = 0; dir < CommonPrms::instance()->Ndim(); ++dir){
      SUNmat Pxm(P_[gf.cslice(0,site, dir)]);
      H_local -= ReTr(Pxm*Pxm);
    }
  }
  double H = Communicator::instance()->reduce_sum(H_local);
  CCIO::cout << "[Momenta] H_p = "<< H << std::endl;

  // action terms
  for(int lv = 0; lv< as_.size(); ++lv)
    for(int id = 0; id < as_.at(lv).size(); ++id)
      H+= as_[lv].at(id)->calc_H();

  return H;
}

void MDexec_leapfrog::
integrator_step(int cl,std::vector<int>& clock){
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
      update_P(cl,eps/2);
      ++clock[cl];
      for(int l=0; l<cl;++l) CCIO::cout<<"   ";
      CCIO::cout<<"P "<< static_cast<double>(clock[cl])/2 <<endl;
    }
    if(cl == fl){          // lowest level 
      update_U(eps);
      for(int l=0; l<cl;++l) CCIO::cout<<"   ";
      CCIO::cout<<"U "<< static_cast<double>(clock[cl]+1)/2 <<endl;
    }else{                 // recursive function call 
      integrator_step(cl+1,clock);
    }
    if(clock[cl] == fin){  // final half step
      update_P(cl,eps/2);
      
      ++clock[cl];
      for(int l=0; l<cl;++l) CCIO::cout<<"   ";
      CCIO::cout<<"P "<< static_cast<double>(clock[cl])/2 <<endl;
    }else{                  // bulk step
      update_P(cl,eps);
      
      clock[cl]+=2;
      for(int l=0; l<cl;++l) CCIO::cout<<"   ";
      CCIO::cout<<"P "<< static_cast<double>(clock[cl])/2 <<endl;
    }
  }
}

void MDexec_leapfrog::
integrator(int cl,std::vector<int>& clock){
  // cl  : current level
  // fl  : final level
  // eps : current step size
  
  for(int step=0; step< Params.MDsteps; ++step){   // MD step 
    //CCIO::cout<<"MDstep = "<< step << endl;
    integrator_step(cl,clock);
  }
}

const Field MDexec_leapfrog::get_U() const{ return *U_;}
