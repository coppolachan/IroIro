/*!
  @file fermion_meas_factory.hpp 
  @brief Fermionic measurements operators factories
*/
#ifndef FERMION_MEAS_ABS_FACT_
#define FERMION_MEAS_ABS_FACT_

#include "include/common_fields.hpp"
#include "quark_propagators.hpp"
#include "include/pugi_interface.h"
/*!
 * @brief Abstract base class for creating QuarkPropagators
 */
class QuarkPropagatorFactory {
public:
  virtual QuarkPropagator* getQuarkProp(GaugeField& Field) = 0;
  virtual ~QuarkPropagatorFactory(){}
};
/*
class GenericMeasurementFactory {
public:
  virtual GenericMeasurement* getGMeas(GaugeField& Field) = 0;
  ~GenericMeasurementFactory();
}
*/

///////////////////////////////////////////////////

namespace QuarkPropagators{
  QuarkPropagatorFactory* createQuarkPropagatorFactory(const XML::node);
}

#endif
