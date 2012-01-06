#ifndef ACTION_GAUGE_FACT_
#define ACTION_GAUGE_FACT_

#include "include/pugi_interface.h"
#include "Action/action_gauge_wilson.hpp"
#include "Action/action_gauge_rect.hpp"
#include "action_Factory.hpp"


class GaugeActionFactory : public ActionFactory {
  virtual Action* getGaugeAction(const Format::Format_G&,
				 Field* const) = 0;
public:
  Action* getAction(const Format::Format_G& GaugeForm,
		    Field* const GaugeField) {
    return getGaugeAction(GaugeForm,
			  GaugeField);
  }
};

///////////////////////////////////////////////////////////////////////

class WilsonGaugeActionFactory : public GaugeActionFactory {
  const XML::node Action_node;
  
public:
  WilsonGaugeActionFactory(const XML::node node):
    Action_node(node){};
  
private:  
  ActionGaugeWilson* getGaugeAction(const Format::Format_G& Form,
				    Field* const GaugeField){
    return new ActionGaugeWilson(Action_node,
				 Form,
				 GaugeField); //pass xml node
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
  ActionGaugeRect* getGaugeAction(const Format::Format_G& Form,
				  Field* const GaugeField){
    return new ActionGaugeRect(Action_node,
			       Form,
			       GaugeField); //pass xml node
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
  ActionGaugeRect* getGaugeAction(const Format::Format_G& Form,
				  Field* const GaugeField){
    return new ActionGaugeRect(Action_node,
			       3.648,
			       -0.331,
			       Form,
			       GaugeField); //pass xml node
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
  ActionGaugeRect* getGaugeAction(const Format::Format_G& Form,
				  Field* const GaugeField){
    return new ActionGaugeRect(Action_node,
			       5.0/3.0,
			       -1.0/12.0,
			       Form,
			       GaugeField); //pass xml node
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
  ActionGaugeRect* getGaugeAction(const Format::Format_G& Form,
				  Field* const GaugeField){
    return new ActionGaugeRect(Action_node,
			       12.2704,
			       -1.4088,
			       Form,
			       GaugeField); //pass xml node
  }
};

//Add new gauge action factories here
//.....

////////////////////////////////////////////////
namespace GaugeAction {
  //we need only one Factory
  static GaugeActionFactory* GaugeAct;
  GaugeActionFactory* createGaugeActionFactory(XML::node);

}

#endif
