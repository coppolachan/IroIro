/*!
  @file action_gauge_wilson.hpp

  @brief Declaration of the ActionGaugeWilson class
*/
#ifndef ACTION_GAUGE_WILSON_INCLUDED
#define ACTION_GAUGE_WILSON_INCLUDED

#include "action.h"
#include "include/common_fields.hpp"
#include "Measurements/GaugeM/staples.h"
#include "include/pugi_interface.h"

/*!
  Parameters for the ActionGaugeWilson class
 */
struct ActionGaugeWilsonPrm {
  double beta;/*!< %Gauge coupling */

  ActionGaugeWilsonPrm(const XML::node node){
    //read parameters from XML tree
    XML::read(node, "beta", beta);
  }
  ActionGaugeWilsonPrm(const double beta_)
    :beta(beta_){}
};

/*!
  Defines a concrete class for %gauge %Actions (Wilson)
 */
class ActionGaugeWilson :public Action {
private:
  ActionGaugeWilsonPrm Params;
  const int Ndim_;
  const int Nvol_;
  const int Nc_;
  const Format::Format_G& gf_;
  Field* const u_;
  const Staples* stpl_;


public:
  void  init(const RandNum&,const void* = 0){
    CCIO::cout<<"ActionGaugeWilson initialized"<<std::endl;
  }
  double calc_H();
  Field  md_force(const void* = 0);

  ActionGaugeWilson(const double beta, 
		    const Format::Format_G& gf, 
		    Field* const GField)
    :u_(GField),
     Params(ActionGaugeWilsonPrm(beta)), 
     Ndim_(CommonPrms::instance()->Ndim()),
     Nvol_(CommonPrms::instance()->Nvol()),
     Nc_(  CommonPrms::instance()->Nc()),
     gf_(gf),
     stpl_(new Staples(gf)){}

  ActionGaugeWilson(const XML::node node, 
		    const Format::Format_G& gf,
		    Field* const GField)
    :u_(GField),
     Params(ActionGaugeWilsonPrm(node)), 
     Ndim_(CommonPrms::instance()->Ndim()),
     Nvol_(CommonPrms::instance()->Nvol()),
     Nc_(  CommonPrms::instance()->Nc()),
     gf_(gf),
     stpl_(new Staples(gf_)){}

  ~ActionGaugeWilson(){}
};
#endif
