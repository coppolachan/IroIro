#ifndef ACTION_GAUGE_FACT_ABS_
#define ACTION_GAUGE_FACT_ABS_

#include "include/pugi_interface.h"
#include "action_Factory.hpp"


class GaugeActionFactory : public ActionFactory {
  virtual Action* getGaugeAction(GaugeField* const) = 0;
public:
  Action* getAction(GaugeField* const G) {
    return getGaugeAction(G);
  }
};
////////////////////////////////////////////////
namespace GaugeAction {
  //we need only one Factory
  static GaugeActionFactory* GaugeAct;
  GaugeActionFactory* createGaugeActionFactory(XML::node);

}

#endif
