/*!
 * @file action_gaugetype_factory_abs.hpp 
 *
 * @brief Declaration of the abstract gauge action factory
 *
 * Time-stamp: <2013-04-22 16:38:28 neo>
 */

#ifndef ACTION_GAUGE_FACT_ABS_
#define ACTION_GAUGE_FACT_ABS_

#include "action_factory.hpp"

/*
 * @brief Abstract class for gauge action creation
 */
class GaugeActionFactory : public ActionFactory {
  virtual Action* getGaugeAction(GaugeField* const, SmartConf* const SC) = 0;
public:
  Action* getAction(GaugeField* const G, SmartConf* const SC) {
    return getGaugeAction(G, SC);
  }
};

#endif
