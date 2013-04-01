/*!
  @file action_gauge_wilson_adjoint.hpp

  @brief Declaration of the ActionGaugeWilsonAdjoint class
*/
#ifndef ACTION_GAUGE_WILSON_ADJ_INCLUDED
#define ACTION_GAUGE_WILSON_ADJ_INCLUDED

#include "action.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include "include/pugi_interface.h"

/*!
  Parameters for the ActionGaugeWilson class
 */
struct ActionGaugeWilsonAdjPrm {
  double beta;/*!< %Gauge coupling */

  ActionGaugeWilsonAdjPrm(const XML::node node){
    XML::read(node, "beta", beta);
  }
  ActionGaugeWilsonAdjPrm(const double beta_)
    :beta(beta_){}
};

/*!
  @brief Defines a concrete class for %gauge %Actions (Wilson)
 */
class ActionGaugeWilsonAdjoint :public Action {
private:
  GaugeField* const u_;
  const Staples stpl_;
  ActionGaugeWilsonAdjPrm Params;
  const int Nvol_;

public:
  void  init(const RandNum&){}
  double calc_H();
  GaugeField  md_force();
  void observer_update(){}
  
  ActionGaugeWilsonAdjoint(const double beta, 
			   GaugeField* const GField)
    :u_(GField),
     stpl_(),
     Params(ActionGaugeWilsonAdjPrm(beta)), 
     Nvol_(CommonPrms::instance()->Nvol()){}
  
  /*!
   * @brief Constructor with XML reader
   *
   */
  ActionGaugeWilsonAdjoint(const XML::node node,
		    GaugeField* const GField)
    :u_(GField),
     stpl_(),
     Params(ActionGaugeWilsonAdjPrm(node)), 
     Nvol_(CommonPrms::instance()->Nvol()){}
  
  ~ActionGaugeWilsonAdjoint(){}
};
#endif
