/*!
  @file action_gauge_rect.hpp

  @brief Declaration of the ActionGaugeRect class
*/
#ifndef ACTION_GAUGE_CLOVER_INCLUDED
#define ACTION_GAUGE_CLOVER_INCLUDED

#include "action.h"
#include "include/common_fields.hpp"
#include "Measurements/GaugeM/staples.h"
#include "include/pugi_interface.h"


/*!
  @brief Parameters for the ActionGauge class

  Iwasaki action: \f$c_1 =  3.648\f$, \f$c_2 = -0.331\f$
  Symanzik action: \f$c_1 =  5/3\f$, \f$c_2 = -1/12\f$
  DBW2 action: \f$c_1 =  12.2704\f$, \f$c_2 = -1.4088\f$
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
  ActionGaugeRectPrm Params;
  const int Ndim_;
  const int Nvol_;
  const int Nc_;
  const Format::Format_G& gf_;
  const Format::Format_G* sf_;
  const Staples* stpl_;
  Field* const u_;

public:
  void  init(const RandNum&,const void* = 0){};
  double calc_H();
  Field  md_force(const void* = 0);
  
  ActionGaugeRect(const double beta, 
		  const double c_plaq_,
		  const double c_rect_, 
		  const Format::Format_G& gf, 
		  Field* const GField)
    :u_(GField),
     Params(ActionGaugeRectPrm(beta, c_plaq_, c_rect_)), 
     Ndim_(CommonPrms::instance()->Ndim()),
     Nvol_(CommonPrms::instance()->Nvol()),
     Nc_(  CommonPrms::instance()->Nc()),
     gf_(gf),
     sf_(new Format::Format_G(Nvol_,1)),
     stpl_(new Staples(gf_)){}
  
  /*!
   * @brief Constructor for factories
   *
   */
  ActionGaugeRect(const XML::node node, 
		  const Format::Format_G& gf,
		  Field* const GField)
    :u_(GField),
     Params(ActionGaugeRectPrm(node)), 
     Ndim_(CommonPrms::instance()->Ndim()),
     Nvol_(CommonPrms::instance()->Nvol()),
     Nc_(  CommonPrms::instance()->Nc()),
     gf_(gf),
     sf_(new Format::Format_G(Nvol_,1)),
     stpl_(new Staples(gf_)){}

  /*!
   * @brief Specialized constructor for the case of
   * prefedined coefficients c_plaq and c_rect
   *
   * To be used in factories
   */
  ActionGaugeRect(const XML::node node, 
		  const double c_plaq,
		  const double c_rect,
		  const Format::Format_G& gf,
		  Field* const GField)
    :u_(GField),
     Params(ActionGaugeRectPrm(node, c_plaq, c_rect)), 
     Ndim_(CommonPrms::instance()->Ndim()),
     Nvol_(CommonPrms::instance()->Nvol()),
     Nc_(  CommonPrms::instance()->Nc()),
     gf_(gf),
     sf_(new Format::Format_G(Nvol_,1)),
     stpl_(new Staples(gf_)){}

  
  ~ActionGaugeRect(){delete stpl_;}
};



#endif
