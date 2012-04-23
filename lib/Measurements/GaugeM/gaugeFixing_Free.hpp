/*! @file gaugeFixing_Free.hpp
    @brief GaugeFixing_Free class defines gauge unfixing
*/

#ifndef GAUGEFIXING_FREE_INCLUDED
#define GAUGEFIXING_FREE_INCLUDED

#include "gaugeFixing.hpp"

class GaugeFixing_Free :public GaugeFixing{
public:
  const GaugeField do_fix(const GaugeField& U)const{ return U; }
};

#endif
