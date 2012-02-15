//----------------------------------------------------------------------
/*! @file staples.hpp
  
  @brief Defines the staples measurement classes

*/
//----------------------------------------------------------------------
#ifndef STAPLES_INCLUDED
#define STAPLES_INCLUDED

#include "include/commonPrms.h"
#include "Communicator/communicator.h"
#include "include/field.h"
#include "include/format_G.h"
#include "Main/Geometry/shiftField.h"
#include "Tools/sunMat.h"


#include "include/common_fields.hpp"
//temporary
#include "lib/Main/Geometry/mapper.hpp"

class Staples {
  int Nc_;
  int Ndim_;
  int Nvol_, Lvol_;
  const Format::Format_G& gf_;
  Mapper shift;

  Communicator* com_;
  SiteIndex* idx_;  
  Format::Format_G* format1d_;


public:
  Staples(const Format::Format_G& gf)
    :gf_(gf),
     Nc_(CommonPrms::instance()->Nc()),
     Ndim_(CommonPrms::instance()->Ndim()),
     Nvol_(CommonPrms::instance()->Nvol()),
     Lvol_(CommonPrms::instance()->Lvol()),
     com_(Communicator::instance()),
     idx_(SiteIndex::instance()),
     shift(),
     format1d_(new Format::Format_G(Nvol_,1)){}

  ~Staples(){delete format1d_;}

  /////////////////////////////////////////////////
  double plaquette(const GaugeFieldType&) const;
  double plaq_s   (const GaugeFieldType&) const;
  double plaq_t   (const GaugeFieldType&) const;
  GaugeField1DType lower(const GaugeFieldType&, int, int) const;

  //////////////////////////////////////////////////

  // functions for Field class
  ShiftField_up<GaugeFieldFormat> LinkUp(const Field&, int, int) const;
  Field upper(const Field&, int, int) const;
  Field lower(const Field&, int, int) const;

  double plaq_s(const Field&) const;
  double plaq_t(const Field&) const;
  double plaquette(const Field&) const;
  
  void staple(Field&, const Field&, int) const;

  // functions for ShiftField class
  Field upper(const ShiftField&, int, int) const;
  Field lower(const ShiftField&, int, int) const;

  double plaq_s(const ShiftField&) const;
  double plaq_t(const ShiftField&) const;
  double plaquette(const ShiftField&) const;
  
  void staple(Field&, const ShiftField&, int) const;
};









#endif  
