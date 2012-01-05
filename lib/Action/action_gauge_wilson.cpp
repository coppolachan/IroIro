/*!
  @file action_gauge_wilson.cpp

  @brief Definition of the ActionGaugeWilson class
*/
#include "action_gauge_wilson.hpp"
#include "Communicator/comm_io.hpp"

using namespace std;


double ActionGaugeWilson::calc_H(){
  using namespace SUNmat_utils;
  double plaq = stpl_->plaquette(*u_);//already reduced
  int NP = CommonPrms::instance()->NP();

  CCIO::cout<<" -- Plaquette = "<<plaq<<endl;
 
  double Hgauge = Params.beta*Nvol_*NP*Ndim_*(Ndim_-1)/2*(1-plaq);
#if VERBOSITY>=ACTION_VERB_LEVEL
  CCIO::cout << "[ActionGaugeWilson] H = "<<Hgauge<<std::endl;
#endif
  return Hgauge;
}


Field ActionGaugeWilson::md_force(const void*){
  using namespace SUNmat_utils;
  SUNmat pl;
  Field force(gf_.size());
  
  for(int m = 0; m < Ndim_; ++m){
    
    GaugeField1D tmp;
    
    for(int n=0; n< Ndim_; ++n){
      if(n != m){
	tmp.U += stpl_->upper(*u_,m,n);
	tmp.U += stpl_->lower(*u_,m,n);
      }
    }
    for(int site=0; site<Nvol_; ++site){
      pl = (u(*u_,gf_,site,m)*u_dag(tmp,site));
      force.set(gf_.cslice(0,site,m), anti_hermite(pl));
    }
  }
  force *= 0.5*Params.beta/Nc_;

#if VERBOSITY>=ACTION_VERB_LEVEL
  monitor_force(force, "ActionGaugeWilson");
#endif

  return force;
}


