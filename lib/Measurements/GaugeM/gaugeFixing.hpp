/*! @file gaugeFixing.hpp
    @brief definition of GaugeFixing interface
*/
#ifndef GAUGEFIXING_INCLUDED
#define GAUGEFIXING_INCLUDED

#include "include/common_fields.hpp"

class GaugeFixing{
public:
  GaugeFixing(){}
  virtual ~GaugeFixing(){}
  virtual const GaugeField do_fix(const GaugeField&)const = 0;
};

#endif
