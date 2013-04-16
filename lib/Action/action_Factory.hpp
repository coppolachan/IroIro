/*!
  @file action_Factory.hpp

  @brief Declaration of abstract ActionFactory class, top level

  Time-stamp: <2013-04-16 11:22:24 neo>
 */


#ifndef ACTION_FACT_HPP_
#define ACTION_FACT_HPP_

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
