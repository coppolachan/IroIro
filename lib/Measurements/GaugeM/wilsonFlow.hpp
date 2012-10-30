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

class WilsonFlow{
private:
  int Nexp_;
  int Nstep_;
  double estep_;
  mutable GaugeField U_; // internal buffer, initialized by input U
  ActionGaugeWilson* Sg_;

  void update_U(const GaugeField& Z) const;
  void flow_step() const;
public:
  WilsonFlow(XML::node node,const GaugeField& U)
    :U_(U),Sg_(new ActionGaugeWilson(3.0,&U_)){
    /*! @brief set beta=3.0 to reuse ActionGaugeWilson */

    XML::read(node,"Nexp", Nexp_, MANDATORY);
    XML::read(node,"Nstep",Nstep_,MANDATORY);
    XML::read(node,"estep",estep_,MANDATORY);
  }

  WilsonFlow(double beta,int Nexp,int Nstep,double estep,const GaugeField& U)
    :Nexp_(Nexp),Nstep_(Nstep),estep_(estep),U_(U),
     Sg_(new ActionGaugeWilson(3.0,&U_)){}
  
  ~WilsonFlow(){if(Sg_) delete Sg_;}
  
  std::vector<double> evolve()const;
  std::vector<double> get_t()const;
};

#endif
