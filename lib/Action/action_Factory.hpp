#ifndef ACTION_FACT_
#define ACTION_FACT_

#include "action.hpp"
#include "Smearing/smearingFactories.hpp"

/*!
 *@class ActionFactory
 *
 *@brief Abstract Factory class for General actions
 *
 */
class ActionFactory {
public:
  /*!
    Virtual function returning the %Action
  */
  virtual Action* getAction(GaugeField* const, SmartConf* const) = 0;
};


#endif
