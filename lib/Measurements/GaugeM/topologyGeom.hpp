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
  void calc_Q(std::vector<double>& q,const GaugeField&)const;
public:
  TopologyGeom():stpl_(),Nvol_(CommonPrms::instance()->Nvol()){}

  double get_Q(const GaugeField&)const;
  double get_Q(std::vector<double>&,const GaugeField&)const;
};

#endif
