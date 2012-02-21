/*!
  @file action_gauge_wilson.hpp

  @brief Declaration of the ActionGaugeWilson class
*/
#ifndef ACTION_GAUGE_WILSON_INCLUDED
#define ACTION_GAUGE_WILSON_INCLUDED

#include "action.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include "include/pugi_interface.h"

/*!
  Parameters for the ActionGaugeWilson class
 */
struct ActionGaugeWilsonPrm {
  double beta;/*!< %Gauge coupling */

  ActionGaugeWilsonPrm(const XML::node node){
    XML::read(node, "beta", beta);
  }
  ActionGaugeWilsonPrm(const double beta_)
    :beta(beta_){}
};

/*!
  @brief Defines a concrete class for %gauge %Actions (Wilson)
 */
class ActionGaugeWilson :public Action {
private:
  GaugeField* const u_;
  const Staples stpl_;
  ActionGaugeWilsonPrm Params;
  const int Ndim_;
  const int Nvol_;

public:
  void  init(const RandNum&,const void* = 0){}
  double calc_H();
  Field  md_force(const void* = 0);
  void observer_update(){}; 

  ActionGaugeWilson(const double beta, 
		    GaugeField* const GField)
    :u_(GField),
     stpl_(),
     Params(ActionGaugeWilsonPrm(beta)), 
     Ndim_(CommonPrms::instance()->Ndim()),
     Nvol_(CommonPrms::instance()->Nvol()){}


  /*!
   * @brief Constructor with XML reader
   *
   */
  ActionGaugeWilson(const XML::node node,
		    GaugeField* const GField)
    :u_(GField),
     stpl_(),
     Params(ActionGaugeWilsonPrm(node)), 
     Ndim_(CommonPrms::instance()->Ndim()),
     Nvol_(CommonPrms::instance()->Nvol()){}

  ~ActionGaugeWilson(){}
};
#endif
