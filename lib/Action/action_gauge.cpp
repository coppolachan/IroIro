/*!
  @file action_gauge.cpp

  @brief Implementation of the Action_gauge class
*/
#include "action_gauge.hpp"
#include "Communicator/comm_io.hpp"

using namespace std;

inline SUNmat ActionGaugeWilson::u(const Field& g,int site,int dir){
  return SUNmat(g[gf_.cslice(0,site,dir)]);
}
inline SUNmat ActionGaugeWilson::u_dag(const Field& g,int site,int dir){
  return SUNmat(g[gf_.cslice(0,site,dir)]).dag();
}

inline SUNmat ActionGaugeWilson::u(const Field& g,int site){
  return SUNmat(g[sf_->cslice(0,site)]);
}
inline SUNmat ActionGaugeWilson::u_dag(const Field& g,int site){
  return SUNmat(g[sf_->cslice(0,site)]).dag();
}

double ActionGaugeWilson::calc_H(){
  //CCIO::cout<<"ActionGaugeWilson::calc_H"<<endl;
  double plaq = stpl_->plaquette(*u_);
  int NP = CommonPrms::instance()->NP();

  CCIO::cout<<"Plaq = "<<plaq<<endl;
 
  double Hgauge = Params.beta*Nvol_*NP*Ndim_*(Ndim_-1)/2*(1-plaq);
  CCIO::cout << "[Action_Gauge] H = "<<Hgauge<<std::endl;
  return Hgauge;
}

Field ActionGaugeWilson::md_force(const void*){
  using namespace SUNmat_utils;
  
  Field force(gf_.size());
  
  for(int m = 0; m < Ndim_; ++m){
    
    Field tmp(sf_->size());
    
    for(int n=0; n< Ndim_; ++n){
      if(n != m){
	tmp += stpl_->upper(*u_,m,n);
	tmp += stpl_->lower(*u_,m,n);
      }
    }
    for(int site=0; site<Nvol_; ++site){
      
      SUNmat staple_sum = u_dag(tmp,site);
      SUNmat link = u(*u_,site,m);
      SUNmat pl = link*staple_sum;
      
      force.set(gf_.cslice(0,site,m), anti_hermite(pl));
    }
  }
  force *= Params.beta/Nc_/2.0;

#if VERBOSITY>=ACTION_VERB_LEVEL
  monitor_force(force, "Action_Gauge");
#endif

  return force;
}


