/*
 * @file mdExec_leapfrog.cpp
 *
 * @brief Definition of MDexec_leapfrog class and Parameters 
 */
#include "mdExec_leapfrog.hpp"
#include "Tools/sunMatUtils.hpp"
#include "include/messages_macros.hpp"

#include <iomanip>

using namespace std;
//using namespace Format;

void MDexec_leapfrog::register_observers(){
  // Register actions 
  for(int level = 0; level < as_.size(); ++level){
    for(int id = 0; id < as_.at(level).size(); ++id){
      _Message(DEBUG_VERB_LEVEL, "Registering Observers - Action level = "
	       << level <<" Action# = "<< id<<"\n");
      attach_observer(GaugeObservers, as_[level].at(id));
    }
  }
  
  // Register other observers
  // .....
  CCIO::cout << "[MDexec_leapfrog] Registered "<<GaugeObservers.size()<<" Gauge observers\n";
}

void MDexec_leapfrog::attach_observer(ObserverList& OList,
				      Observer* Obs){
  OList.push_back(Obs);
}

void MDexec_leapfrog::notify_observers(ObserverList& OList) {
  for(int element = 0; element < OList.size(); ++element) {
    OList[element]->observer_update();
  }
}


void MDexec_leapfrog::update_U(double ep){
  using namespace SUNmatUtils;
  using namespace FieldUtils;

  int Ndim = CommonPrms::instance()->Ndim();
  int Nvol = CommonPrms::instance()->Nvol();

  const SUNmat I = unity();
  SUNmat au;

  for(int m = 0; m < Ndim; ++m){
    for(int site=0; site<Nvol; ++site){
      au = I;
      for(int k = Params.Nexp; k > 0; --k){
 	au *= ep/k;
	au *= matrix(P_,site,m);
	au += I;
      }
      au *= matrix(*U_,site,m);
      U_->data.set(U_->format.islice(site,m),au.reunit().getva());
    }
  }

  notify_observers(GaugeObservers);
}

void MDexec_leapfrog::update_P(int lv,double ep){
  for(int a = 0; a < as_[lv].size(); ++a){
    GaugeField fce = as_[lv].at(a)->md_force();
    fce *= ep;
    P_-= fce; 
  }
}

void MDexec_leapfrog::
init(vector<int>& clock,const GaugeField& U,const RandNum& rand){
  clock.resize(as_.size(),0.0);  

  *U_= U;                       // initialize U_ (common to actions) to U
  notify_observers(GaugeObservers);
  MDutils::md_mom(P_,rand); // initialize P_ 

  for(int lv = 0; lv< as_.size(); ++lv){
    for(int id = 0; id < as_.at(lv).size(); ++id){
      _Message(DEBUG_VERB_LEVEL, "Initialization of MD steps level = "<< 
	       lv <<" Action# = "<< id<<"\n");
      as_[lv].at(id)->init(rand);
    }
  }



}

double MDexec_leapfrog::calc_H()const{
  using namespace SUNmatUtils;
  // kinetic term
  double H_local = 0.0;

  for(int site = 0; site < CommonPrms::instance()->Nvol(); ++site){
    for(int dir = 0; dir < CommonPrms::instance()->Ndim(); ++dir){
      SUNmat Pxm(P_.data[P_.format.cslice(0,site, dir)]);
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
    integrator_step(cl,clock);
  }
}

const GaugeField MDexec_leapfrog::get_U() const{ return *U_;}
