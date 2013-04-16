/*!
  @file action_gauge_wilson_force.cpp
  @brief Definition of the md_force method for the ActionGaugeWilson class

  This is the standard, plain c++ version

  Time-stamp: <2013-04-16 16:17:02 neo>
*/
#include "action_gauge_wilson.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/fieldUtils.hpp"
#include "include/messages_macros.hpp"

GaugeField ActionGaugeWilson::md_force(){
  using namespace FieldUtils;
  using namespace SUNmatUtils;

  SUNmat pl;
  GaugeField force;
  GaugeField1D tmp; 
 
  for(int m = 0; m < NDIM_; ++m){
    tmp = 0.0;
    for(int n=0; n< NDIM_; ++n)
      if(n != m) tmp += stpl_.upper_lower(*u_,m,n);

    for(int site=0; site < Nvol_; ++site){
      pl = mat(*u_,site,m)*mat_dag(tmp,site);
      SetMat(force, anti_hermite_traceless(pl), site, m);
    }
  }

  force *= 0.5*Params.beta/NC_;

  _MonitorMsg(ACTION_VERB_LEVEL, Action, force, "ActionGaugeWilson");

  return force;
}


