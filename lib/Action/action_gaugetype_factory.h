#ifndef ACTION_GAUGE_FACT_
#define ACTION_GAUGE_FACT_

#include "include/pugi_interface.h"
#include "Action/action_gauge.h"
#include "action_Factory.h"


class GaugeActionCreator : public ActionCreator {
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

class WilsonGaugeActionCreator : public GaugeActionCreator {
  const XML::node Action_node;
  
public:
  WilsonGaugeActionCreator(const XML::node node):
    Action_node(node){};
  
private:  
  ActionGaugeWilson* getGaugeAction(const Format::Format_G& Form,
				    Field* const GaugeField){
    return new ActionGaugeWilson(Action_node,
				 Form,
				 GaugeField); //pass xml node
  };
  
};

//Add new gauge action factories here
//.....

////////////////////////////////////////////////
namespace GaugeAction {
  //we need only one Creator
  static GaugeActionCreator* GaugeAct;
  GaugeActionCreator* createGaugeActionFactory(XML::node);

}

#endif
