/*! @file gaugeFixing_Free.hpp
    @brief GaugeFixing_Free class defines gauge unfixing
*/

#ifndef GAUGEFIXING_FREE_INCLUDED
#define GAUGEFIXING_FREE_INCLUDED

#include "include/common_fields.hpp"

class GaugeFixing_Free :public GaugeFixing{
public:
  const GaugeField fix(const GaugeField& U)const{ return U; }
};

#endif
