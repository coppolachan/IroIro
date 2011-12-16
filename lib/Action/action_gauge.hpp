/*!
  @file action_gauge.h

  @brief Declaration of the Action_gauge class
*/
#ifndef ACTION_GAUGE_INCLUDED
#define ACTION_GAUGE_INCLUDED

#include "include/common_fields.hpp"
#include "action.h"
#include "Measurements/GaugeM/staples.h"
#include "include/pugi_interface.h"

/*!
  Parameters for the ActionGauge class
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
  Defines a concrete class for %gauge %Actions
 */
class ActionGaugeWilson :public Action {
private:
  ActionGaugeWilsonPrm Params;
  const int Ndim_;
  const int Nvol_;
  const int Nc_;
  const Format::Format_G& gf_;
  const Format::Format_G* sf_;
  const Staples* stpl_;
  Field* const u_;
  int nodeid_;

  SUNmat u(const Field&, int site, int dir);
  SUNmat u(const Field&, int site);

  SUNmat u_dag(const Field&, int site, int dir);
  SUNmat u_dag(const Field&, int site);

public:
  void  init(const RandNum&,const void* = 0){
    CCIO::cout<<"Action_gauge initialized"<<std::endl;
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
     sf_(new Format::Format_G(Nvol_,1)),
     stpl_(new Staples(gf_)),
     nodeid_(Communicator::instance()->nodeid()){}

  ActionGaugeWilson(const XML::node node, 
		    const Format::Format_G& gf,
		    Field* const GField)
    :u_(GField),
     Params(ActionGaugeWilsonPrm(node)), 
     Ndim_(CommonPrms::instance()->Ndim()),
     Nvol_(CommonPrms::instance()->Nvol()),
     Nc_(  CommonPrms::instance()->Nc()),
     gf_(gf),
     sf_(new Format::Format_G(Nvol_,1)),
     stpl_(new Staples(gf_)),
     nodeid_(Communicator::instance()->nodeid()){}

  ~ActionGaugeWilson(){}
};
#endif
