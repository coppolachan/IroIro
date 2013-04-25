/*! @file staples.hpp
  @brief Defines the staples measurement classes
*/
#ifndef STAPLES_INCLUDED
#define STAPLES_INCLUDED

#include "include/commonPrms.h"
#include "include/common_fields.hpp"
#include "Communicator/communicator.hpp"
#include "Geometry/shiftField.hpp"

class Staples {
  int Nvol_, Lvol_;
  Communicator* com_;

public:
  Staples():Nvol_(CommonPrms::instance()->Nvol()),
	    Lvol_(CommonPrms::instance()->Lvol()),
	    com_(Communicator::instance()){
    Mapping::init_shiftField();
  }

  double plaquette(const GaugeField&) const;
  double plaq_s   (const GaugeField&) const;
  double plaq_t   (const GaugeField&) const;
  double plaq_min (const GaugeField&,double threshold=1.0) const;

  double plaquette_adj(const GaugeField&) const;
  double plaq_s_adj   (const GaugeField&) const;
  double plaq_t_adj   (const GaugeField&) const;

  GaugeField1D lower(const GaugeField&,int,int) const;
  GaugeField1D upper(const GaugeField&,int,int) const;
  GaugeField1D upper_lower(const GaugeField&,int,int) const;
  GaugeField1D upper_lower(const GaugeField&,int,int, const Field) const;
  GaugeField1D fieldStrength(const GaugeField&,int,int)const;

  void staple(GaugeField1D&, const GaugeField&,int) const;
};

#endif  
