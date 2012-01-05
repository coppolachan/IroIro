//----------------------------------------------------------------------
/*!
  @file rectangular.hpp

  @brief Declares RectangularStaple class

*/ 
//----------------------------------------------------------------------
#ifndef RECTANGULAR_INCLUDED
#define RECTANGULAR_INCLUDED

#ifndef COMMONPRMS_INCLUDED
#include "include/commonPrms.h"
#endif

#ifndef COMMUNICATOR_INCLUDED
#include "Communicator/communicator.h"
#endif

#ifndef FIELD_INCLUDED
#include "include/field.h"
#endif

#ifndef FORMAT_G_INCLUDED
#include "include/format_G.h"
#endif

#ifndef SHIFTFIELD_INCLUDED
#include "Main/Geometry/shiftField.h"
#endif

#ifndef SUNMAT_INCLUDED
#include "Tools/sunMat.h"
#endif

#ifndef STAPLES_INCLUDED
#include "staples.h"
#endif

class RectangularStaples: {
private:
  int Nc_;
  int Ndim_;
  int Nvol_, Lvol_;
  const Format::Format_G& gf_;

  Communicator* com_;
  SiteIndex* idx_;  
  Format::Format_G* format1d_;

public:
  RectangularStaples(const Format::Format_G& gf)
    :gf_(gf),
     Nc_(CommonPrms::instance()->Nc()),
     Ndim_(CommonPrms::instance()->Ndim()),
     Nvol_(CommonPrms::instance()->Nvol()),
     Lvol_(CommonPrms::instance()->Lvol()),
     com_(Communicator::instance()),
     idx_(SiteIndex::instance()),
     format1d_(new Format::Format_G(Nvol_,1)){}

  ~Staples(){delete format1d_;}

  // functions for Field class
  Field upper_1(const Field&, int, int) const;
  Field lower_1(const Field&, int, int) const;

  Field upper_2(const Field&, int, int) const;
  Field lower_2(const Field&, int, int) const;

  Field upper_3(const Field&, int, int) const;
  Field lower_3(const Field&, int, int) const;

  // functions for ShiftField class
  Field upper_1(const ShiftField&, int, int) const;
  Field lower_1(const ShiftField&, int, int) const;

  Field upper_2(const ShiftField&, int, int) const;
  Field lower_2(const ShiftField&, int, int) const;

  Field upper_3(const ShiftField&, int, int) const;
  Field lower_3(const ShiftField&, int, int) const;

};

#endif  
