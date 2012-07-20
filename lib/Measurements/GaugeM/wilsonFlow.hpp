/*!@file wilsonFlow.hpp
 * @brief definition of the WilsonFlow class
 */
#ifndef WILSONFLOW_INCLUDED
#define WILSONFLOW_INCLUDED

#include "include/messages_macros.hpp"
#include "include/commonPrms.h"
#include "include/common_fields.hpp"
#include "include/pugi_interface.h"
#include "Action/action_gauge_wilson.hpp"

class WilsonFlow{
private:
  int Lvol_,Nvol_;
  int Nexp_;
  int Nstep_;
  double estep_;
  mutable GaugeField U_; // internal buffer, initialized by input U
  ActionGaugeWilson* Sg_;

  void update_U(const GaugeField& Z) const;
  void gradFlow() const;
public:
  WilsonFlow(XML::node node,const GaugeField& U)
    :U_(U),
     Sg_(new ActionGaugeWilson(NC_,&U_)),
     Nvol_(CommonPrms::instance()->Nvol()),
     Lvol_(CommonPrms::instance()->Lvol()){
    XML::read(node,"Nexp",Nexp_,MANDATORY);
    XML::read(node,"Nstep",Nstep_,MANDATORY);
    XML::read(node,"estep",estep_,MANDATORY);
  }

  WilsonFlow(int Nexp,int Nstep,double estep,const GaugeField& U)
    :Nexp_(Nexp),Nstep_(Nstep),estep_(estep),U_(U),
     Sg_(new ActionGaugeWilson(NC_,&U_)),
     Nvol_(CommonPrms::instance()->Nvol()),
     Lvol_(CommonPrms::instance()->Lvol()){}
  
  ~WilsonFlow(){if(Sg_) delete Sg_;}

  std::vector<double> evolve()const;
};

#endif
