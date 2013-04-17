/*!
 * @file action_gaugetype_factory_creator.hpp 
 *
 * @brief Definition of the GaugeAction factory namespace
 *
 * Time-stamp: <2013-04-17 15:00:06 neo>
 */

#ifndef ACTION_GAUGE_FACT_CREATOR_
#define ACTION_GAUGE_FACT_CREATOR_

#include "include/pugi_interface.h"
#include "action_gaugetype_factory_abs.hpp"

namespace GaugeAction {
  GaugeActionFactory* createGaugeActionFactory(XML::node);
}

#endif
