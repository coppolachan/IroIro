/*! @file staples.hpp
  
  @brief Defines the staples measurement classes
*/
#ifndef STAPLES_INCLUDED
#define STAPLES_INCLUDED

#include "include/commonPrms.h"
#include "include/common_fields.hpp"
#include "Communicator/communicator.h"

class Staples {
  int Nvol_, Lvol_;

  Communicator* com_;

public:
  Staples()
    :Nvol_(CommonPrms::instance()->Nvol()),
     Lvol_(CommonPrms::instance()->Lvol()),
     com_(Communicator::instance()){}

  /////////////////////////////////////////////////
  double plaquette(const GaugeField&) const;
  double plaq_s   (const GaugeField&) const;
  double plaq_t   (const GaugeField&) const;
  GaugeField1D lower(const GaugeField&, int, int) const;
  GaugeField1D upper(const GaugeField&, int, int) const;
  void staple(GaugeField1D&, const GaugeField&, int) const;
  //////////////////////////////////////////////////

};

#endif  
