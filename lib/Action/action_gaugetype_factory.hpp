#ifndef ACTION_GAUGE_FACT_
#define ACTION_GAUGE_FACT_

#include "include/pugi_interface.h"
#include "action_gaugetype_factory_abs.hpp"
#include "Action/action_gauge_wilson.hpp"
#include "Action/action_gauge_wilson_adjoint.hpp"
#include "Action/action_gauge_rect.hpp"


///////////////////////////////////////////////////////////////////////
class WilsonGaugeActionFactory : public GaugeActionFactory {
  const XML::node Action_node;
  
public:
  WilsonGaugeActionFactory(const XML::node node):
    Action_node(node){};
  
private:  
  ActionGaugeWilson* getGaugeAction(GaugeField* const G, SmartConf* const SC){
    return new ActionGaugeWilson(Action_node,G); 
  }
};

///////////////////////////////////////////////////////////////////////
class WilsonGaugeAdjointActionFactory : public GaugeActionFactory {
  const XML::node Action_node;
  
public:
  WilsonGaugeAdjointActionFactory(const XML::node node):
    Action_node(node){};
  
private:  
  ActionGaugeWilsonAdjoint* getGaugeAction(GaugeField* const G, SmartConf* const SC){
    return new ActionGaugeWilsonAdjoint(Action_node,G); 
  }
};

///////////////////////////////////////////////////////////////////////
/*!
 * @ brief Generic rectangle action
 *
 */
class RectGaugeActionFactory : public GaugeActionFactory {
  const XML::node Action_node;
  
public:
  RectGaugeActionFactory(const XML::node node):
    Action_node(node){};
  
private:  
  ActionGaugeRect* getGaugeAction(GaugeField* const G, SmartConf* const SC){
    return new ActionGaugeRect(Action_node,G); 
  }
};
///////////////////////////////////////////////////////////////////////
/*!
 * @ brief Iwasaki action \f$c_plaq =  3.648\f$, \f$c_rect = -0.331\f$
 *
 */
class IwasakiGaugeActionFactory : public GaugeActionFactory {
  const XML::node Action_node;
  
public:
  IwasakiGaugeActionFactory(const XML::node node):
    Action_node(node){};
  
private:  
  ActionGaugeRect* getGaugeAction(GaugeField* const G, SmartConf* const SC){
    return new ActionGaugeRect(Action_node,
			       3.648,
			       -0.331,
			       G);
  }
};
///////////////////////////////////////////////////////////////////////
/*!
 * @ brief Symanzik action \f$c_plaq =  5/3\f$, \f$c_rect = -1/12\f$
 *
 */
class SymanzikGaugeActionFactory : public GaugeActionFactory {
  const XML::node Action_node;
  
public:
  SymanzikGaugeActionFactory(const XML::node node):
    Action_node(node){};
  
private:  
  ActionGaugeRect* getGaugeAction(GaugeField* const G, SmartConf* const SC){
    return new ActionGaugeRect(Action_node,
			       5.0/3.0,
			       -1.0/12.0,
			       G); //pass xml node
  }
};
///////////////////////////////////////////////////////////////////////
/*!
 * @ brief DBW2 action  \f$c_plaq =  12.2704\f$, \f$c_rect = -1.4088\f$
 *
 */
class DBW2GaugeActionFactory : public GaugeActionFactory {
  const XML::node Action_node;
  
public:
  DBW2GaugeActionFactory(const XML::node node):
    Action_node(node){};
  
private:  
  ActionGaugeRect* getGaugeAction(GaugeField* const G, SmartConf* const SC){
    return new ActionGaugeRect(Action_node,
			       12.2704,
			       -1.4088,
			       G); //pass xml node
  }
};

//Add new gauge action factories here
//.....

#endif
