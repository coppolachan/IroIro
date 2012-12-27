/*! @file topologyGeom.hpp
 *  @brief definition of the TopologyGeom class
 */
#ifndef TOPOLOGYGEOM_INCLUDED
#define TOPOLOGYGEOM_INCLUDED

#include "staples.hpp"

class TopologyGeom {
private:
  Staples stpl_;
  int Nvol_;
public:
  TopologyGeom():stpl_(),Nvol_(CommonPrms::instance()->Nvol()){}
  double getQ(const GaugeField&)const;
};

#endif
