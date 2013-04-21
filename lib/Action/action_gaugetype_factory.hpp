/*!
 * @file action_gaugetype_factory.hpp 
 *
 * @brief Declaration of Gauge-type action factories
 *
 * Time-stamp: <2013-04-17 15:13:49 neo>
 */

#ifndef ACTION_GAUGE_FACT_
#define ACTION_GAUGE_FACT_


#include "include/pugi_interface.h"

// Include here the files containing 
// the declarations of the actions
#include "action_gaugetype_factory_abs.hpp"
#include "Action/action_gauge_wilson.hpp"
#include "Action/action_gauge_wilson_adjoint.hpp"
#include "Action/action_gauge_rect.hpp"


///////////////////////////////////////////////////////////////////////
/*!
 * @brief Standard Wilson action
 *
 */
class WilsonGaugeActionFactory : public GaugeActionFactory {
  const XML::node Action_node;

  ActionGaugeWilson* getGaugeAction(GaugeField* const,
				    SmartConf* const);  
public:
  WilsonGaugeActionFactory(const XML::node);
};

///////////////////////////////////////////////////////////////////////
/*!
 * @brief Wilson adjoint links action
 *
 */
class WilsonGaugeAdjointActionFactory : public GaugeActionFactory {
  const XML::node Action_node;

  ActionGaugeWilsonAdjoint* getGaugeAction(GaugeField* const,
					   SmartConf* const);  
public:
  WilsonGaugeAdjointActionFactory(const XML::node node);
};

///////////////////////////////////////////////////////////////////////
/*!
 * @brief Generic rectangle action
 *
 */
class RectGaugeActionFactory : public GaugeActionFactory {
  const XML::node Action_node;

  ActionGaugeRect* getGaugeAction(GaugeField* const,
				  SmartConf* const);
public:
  RectGaugeActionFactory(const XML::node);
};
///////////////////////////////////////////////////////////////////////
/*!
 * @brief Iwasaki action \f$c_plaq =  3.648\f$, \f$c_rect = -0.331\f$
 *
 */
class IwasakiGaugeActionFactory : public GaugeActionFactory {
  const XML::node Action_node;
  
  ActionGaugeRect* getGaugeAction(GaugeField* const,
				  SmartConf* const);  
public:
  IwasakiGaugeActionFactory(const XML::node);  
  
};
///////////////////////////////////////////////////////////////////////
/*!
 * @brief Symanzik action \f$c_plaq =  5/3\f$, \f$c_rect = -1/12\f$
 *
 */
class SymanzikGaugeActionFactory : public GaugeActionFactory {
  const XML::node Action_node;

  ActionGaugeRect* getGaugeAction(GaugeField* const, 
				  SmartConf* const);
public:
  SymanzikGaugeActionFactory(const XML::node node);  
};
///////////////////////////////////////////////////////////////////////
/*!
 * @brief DBW2 action  \f$c_plaq =  12.2704\f$, \f$c_rect = -1.4088\f$
 *
 */
class DBW2GaugeActionFactory : public GaugeActionFactory {
  const XML::node Action_node;

  ActionGaugeRect* getGaugeAction(GaugeField* const, 
				  SmartConf* const);
public:
  DBW2GaugeActionFactory(const XML::node node);  
};

//Add new gauge action factories here
//.....

#endif
