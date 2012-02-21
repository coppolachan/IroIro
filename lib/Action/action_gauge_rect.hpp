/*!
  @file action_gauge_rect.hpp

  @brief Declaration of the ActionGaugeRect class
*/
#ifndef ACTION_GAUGE_CLOVER_INCLUDED
#define ACTION_GAUGE_CLOVER_INCLUDED

#include "action.hpp"
#include "include/common_fields.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include "include/pugi_interface.h"


/*!
  @brief Parameters for the ActionGauge class

  - Iwasaki action: \f$c_1 =  3.648\f$, \f$c_2 = -0.331\f$
  - Symanzik action: \f$c_1 =  5/3\f$, \f$c_2 = -1/12\f$
  - DBW2 action: \f$c_1 =  12.2704\f$, \f$c_2 = -1.4088\f$
 */
struct ActionGaugeRectPrm {
  double beta;/*!< %Gauge coupling */
  double c_plaq;/*!< Coefficient for the plaquette term */
  double c_rect;/*!< Coefficient for the rectangular term */

  ActionGaugeRectPrm(const XML::node node){
    //read parameters from XML tree
    XML::read(node, "beta", beta, MANDATORY);
    XML::read(node, "c_plaq", c_plaq, MANDATORY);
    XML::read(node, "c_rect", c_rect, MANDATORY);
  }

  /*! 
   * @brief Constructor to accomodate the case 
   * of predefined coefficients in factories 
   *
   */
  ActionGaugeRectPrm(const XML::node node,
		     const double c_plaq_,
		     const double c_rect_)
    :c_plaq(c_plaq_),c_rect(c_rect_)
  {
    XML::read(node, "beta", beta, MANDATORY);
  }
  
  ActionGaugeRectPrm(const double beta_,
		     const double c_plaq_,
		     const double c_rect_)
    :beta(beta_),
     c_plaq(c_plaq_),
     c_rect(c_rect_){}
};

/*!
  @brief Defines a concrete class for %gauge %Actions (Rect Type)
  
  Gauge action with plaquette and rectangular Wilson loops.
  Iwasaki, Luscher-Weisz, DBW2 are examples of this type
  of action.
 */
class ActionGaugeRect :public Action {
private:
  GaugeField* const u_;
  const Staples stpl_;
  ActionGaugeRectPrm Params;
  const int Nvol_;

public:
  void  init(const RandNum&){};
  double calc_H();
  GaugeField  md_force();
  void observer_update(){};
  
  ActionGaugeRect(const double beta, 
		  const double c_plaq_,
		  const double c_rect_,  
		  GaugeField* const GField)
    :u_(GField),
     stpl_(),
     Params(ActionGaugeRectPrm(beta, c_plaq_, c_rect_)), 
     Nvol_(CommonPrms::instance()->Nvol()){}
  
  /*!
   * @brief Constructor for factories
   *
   */
  ActionGaugeRect(const XML::node node,
		  GaugeField* const GField)
    :u_(GField),
     stpl_(),
     Params(ActionGaugeRectPrm(node)), 
     Nvol_(CommonPrms::instance()->Nvol()){}

  /*!
   * @brief Specialized constructor for the case of
   * prefedined coefficients c_plaq and c_rect
   *
   * To be used in factories
   */
  ActionGaugeRect(const XML::node node, 
		  const double c_plaq,
		  const double c_rect,
		  GaugeField* const GField)
    :u_(GField),
     stpl_(),
     Params(ActionGaugeRectPrm(node, c_plaq, c_rect)), 
     Nvol_(CommonPrms::instance()->Nvol()){}
  
};
#endif
