#ifndef ACTION_FACT_
#define ACTION_FACT_

#include "action.hpp"
#include "include/format_G.h"

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
  virtual Action* getAction(GaugeField* const) = 0;
};


#endif
