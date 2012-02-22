/*
  @file BaseSmear.hpp

  @brief Declares base smearing class Smear

 */

#ifndef BASE_SMEAR_H
#define BASE_SMEAR_H

#include "include/common_fields.hpp"


class Smear {
public:
  virtual ~Smear(){}
  
  virtual void smear     (GaugeField&, const GaugeField&) const = 0;
  virtual void derivative(GaugeField&, const GaugeField&, const GaugeField&) const = 0;
};


#endif
