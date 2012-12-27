/*!@file wilsonFlow.hpp
 * @brief definition of the WilsonFlow class
 */
#ifndef WILSONFLOW_INCLUDED
#define WILSONFLOW_INCLUDED

#include "include/messages_macros.hpp"
#include "include/commonPrms.h"
#include "include/common_fields.hpp"
#include "include/pugi_interface.h"
#include "include/macros.hpp"
#include "Action/action_gauge_wilson.hpp"
#include <string>

class WilsonFlow{
private:
  mutable GaugeField U_; // internal buffer, initialized by input U
  int Nexp_;
  int Nstep_;
  double estep_;
  bool saveConf_;

  ActionGaugeWilson* Sg_; 
  void update_U(const GaugeField& Z) const;

public:
  WilsonFlow(XML::node node,const GaugeField& U)
    :U_(U),Sg_(new ActionGaugeWilson(3.0,&U_)),saveConf_(false){
    /*! @brief set beta=3.0 to reuse ActionGaugeWilson */
    if(!XML::attribute_compare(node,"SaveConf","True")) saveConf_=true;

    XML::read(node,"Nexp", Nexp_, MANDATORY);
    XML::read(node,"Nstep",Nstep_,MANDATORY);
    XML::read(node,"estep",estep_,MANDATORY);
  }

  WilsonFlow(double beta,int Nexp,int Nstep,double estep,
	     const GaugeField& U,bool saveConf=false)
    :Nexp_(Nexp),Nstep_(Nstep),estep_(estep),U_(U),
     Sg_(new ActionGaugeWilson(3.0,&U_)),saveConf_(saveConf){}
  
  ~WilsonFlow(){if(Sg_) delete Sg_;}

  const GaugeField& getU()const{ return U_;}
  
  int Nstep(){return Nstep_;}
  
  double tau(int t)const{  return estep_*(t+1.0);}
  double Edens_plaq(int t)const;
  double Edens_clover(int t)const;
  void evolve_step()const;
  void save_config(const std::string&)const;
};

#endif
