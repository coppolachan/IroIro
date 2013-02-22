/*
 * @file mdExec_2MN.cpp
 * @brief Definition of MDexec_2MN class and Parameters 
 */
#include "mdExec_2MN.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/fieldUtils.hpp"
#include "include/messages_macros.hpp"

#include <iomanip>

#ifdef IBM_BGQ_WILSON
//#include <omp.h>
//#include "bgqthread.h"
#include "include/bgq_su3algebra.h"
#endif

using namespace std;

void MDexec_2MN::register_observers(){
  // Register actions 
  for(int level=0; level<as_.size(); ++level){
    for(int id=0; id<as_.at(level).size(); ++id){
      _Message(DEBUG_VERB_LEVEL, "Registering Observers - Action level = "
	       << level <<" Action# = "<< id<<"\n");
      observers_.push_back(as_[level].at(id));
    }
  }
  // Register other observers
  // .....
  CCIO::cout << "[MDexec_2MN] Registered "
	     <<observers_.size()<<" Gauge observers\n";
}

void MDexec_2MN::notify_observers(){
  for(int elem=0; elem<observers_.size(); ++elem)
    observers_[elem]->observer_update();
}

void MDexec_2MN::update_U(double ep){
  using namespace FieldUtils;
#ifdef IBM_BGQ_WILSON
  Exponentiate_BGQ(*U_, P_, ep, Params.Nexp);
#else    
  using namespace SUNmatUtils;
  int Ndim = CommonPrms::instance()->Ndim();
  int Nvol = CommonPrms::instance()->Nvol();

  for(int m=0; m<Ndim; ++m){
    for(int site=0; site<Nvol; ++site){
      SUNmat au = exponential(mat(P_,site,m)*ep,Params.Nexp);
      au *= mat(*U_,site,m);
      U_->data.set(U_->format.islice(site,m),au.reunit().getva());
    }
  }
#endif  
  notify_observers();
}

void MDexec_2MN::update_P(int lv,double ep){
  for(int a=0; a<as_[lv].size(); ++a){
    GaugeField fce = as_[lv].at(a)->md_force();
#ifdef IBM_BGQ_WILSON
    BGWilsonLA_MatMultScalar((__complex__ double*)fce.data.getaddr(0),
			     -ep,
			     CommonPrms::instance()->Nvol()*CommonPrms::instance()->Ndim());
    BGWilsonMatLA_Add((__complex__ double*)P_.data.getaddr(0), 
		      (__complex__ double*)fce.data.getaddr(0),
		      CommonPrms::instance()->Nvol()*CommonPrms::instance()->Ndim()); 
#else
    fce *= ep;
    P_-= fce; 
#endif
  }
}

void MDexec_2MN::
init(vector<int>& clock,const GaugeField& U,const RandNum& rand){
  clock.resize(as_.size(),0.0);  

  *U_= U;                       // initialize U_ (common to actions) to U
  notify_observers();
  MDutils::md_mom(P_,rand); // initialize P_ 

  for(int lv=0; lv< as_.size(); ++lv){
    for(int id=0; id<as_.at(lv).size(); ++id){
      _Message(DEBUG_VERB_LEVEL, "Initialization of MD steps level = "
	       << lv <<" Action# = "<< id<<"\n");
      as_[lv].at(id)->init(rand);
    }
  }
}

double MDexec_2MN::calc_H()const{
  using namespace SUNmatUtils;
  using namespace FieldUtils;
  // kinetic term
  double H_local = 0.0;
  int Nvol = CommonPrms::instance()->Nvol();
  int Ndim = CommonPrms::instance()->Ndim();
  
#ifdef IBM_BGQ_WILSON
  double *P_ptr = const_cast< Field* >(&(P_.data))->getaddr(0);
  for(int site=0; site<Nvol; ++site){
      for(int dir=0; dir< Ndim; ++dir){
	
	int ref = 18*(site+Nvol*dir);
	H_local -= P_ptr[   ref]*P_ptr[    ref];
	H_local += P_ptr[ 1+ref]*P_ptr[ 1+ ref];
	
	H_local -= P_ptr[ 8+ ref]*P_ptr[ 8+ ref];
	H_local += P_ptr[ 9+ ref]*P_ptr[ 9+ ref];
	
	H_local -= P_ptr[16+ ref]*P_ptr[16+ ref];
	H_local += P_ptr[17+ ref]*P_ptr[17+ ref];
	
	H_local -= 2*(P_ptr[ 2+ref]*P_ptr[ 6+ ref]-
		      P_ptr[ 3+ref]*P_ptr[ 7+ ref]+
		      P_ptr[ 4+ref]*P_ptr[12+ ref]-
		      P_ptr[ 5+ref]*P_ptr[13+ ref]+
		      P_ptr[10+ref]*P_ptr[14+ ref]-
		      P_ptr[11+ref]*P_ptr[15+ ref]);
	
      }
    }
 
#else
  for(int site=0; site<CommonPrms::instance()->Nvol(); ++site){
    for(int dir=0; dir<CommonPrms::instance()->Ndim(); ++dir){
      SUNmat Pxm(P_.data[P_.format.cslice(0,site, dir)]);
      H_local -= ReTr(Pxm*Pxm);
    }
  }
#endif
  double H = Communicator::instance()->reduce_sum(H_local);
  CCIO::cout << "[Momenta] H_p = "<< H << std::endl;

  // action terms
  for(int lv=0; lv<as_.size(); ++lv)
    for(int id=0; id<as_.at(lv).size(); ++id)
      H+= as_[lv].at(id)->calc_H();

  return H;
}

void MDexec_2MN::
integrator_step(int cl,std::vector<int>& clock){
  // cl  : current level
  // fl  : final level
  // eps : current step size

  int fl = as_.size() -1;
  double eps = Params.step_size;
  
  for(int l=0; l<=cl; ++l) eps/= 2*Nrel_[l];
  
  int fin = Nrel_[0];
  for(int l=1; l<=cl; ++l) fin*= 2*Nrel_[l];
  fin = 3*Params.MDsteps*fin -1;
  

  for(int e=0; e<Nrel_[cl]; ++e){
    
    if(clock[cl] == 0){    // initial half step 
      update_P(cl,Params.lambda*eps);
      ++clock[cl];
      for(int l=0; l<cl;++l) CCIO::cout<<"   ";
      CCIO::cout<<"P "<< clock[cl] <<endl;
    }

    if(cl == fl){          // lowest level 
      update_U(0.5*eps);

      for(int l=0; l<cl;++l) CCIO::cout<<"   ";
      CCIO::cout<<"U "<< (clock[cl]+1) <<endl;
    }else{                 // recursive function call 
      integrator_step(cl+1,clock);
    }

    update_P(cl,(1.0-2.0*Params.lambda)*eps);
    ++clock[cl];
    for(int l=0; l<cl;++l) CCIO::cout<<"   ";
    CCIO::cout<<"P "<< (clock[cl]) <<endl;

    if(cl == fl){          // lowest level 
      update_U(0.5*eps);
      
      for(int l=0; l<cl;++l) CCIO::cout<<"   ";
      CCIO::cout<<"U "<< (clock[cl]+1) <<endl;
    }else{                 // recursive function call 
      integrator_step(cl+1,clock);
    }    
    



    if(clock[cl] == fin){  // final half step
      update_P(cl,Params.lambda*eps);
      
      ++clock[cl];
      for(int l=0; l<cl;++l) CCIO::cout<<"   ";
      CCIO::cout<<"P "<< clock[cl] <<endl;
    }else{                  // bulk step
      update_P(cl,Params.lambda*2.0*eps);
      
      clock[cl]+=2;
      for(int l=0; l<cl;++l) CCIO::cout<<"   ";
      CCIO::cout<<"P "<< clock[cl] <<endl;
    }
  }
}

void MDexec_2MN::
integrator(int cl,std::vector<int>& clock){
  // cl  : current level
  // fl  : final level
  // eps : current step size
  
  for(int step=0; step< Params.MDsteps; ++step)   // MD step 
    integrator_step(cl,clock);
}

const GaugeField MDexec_2MN::get_U() const{ return *U_;}
