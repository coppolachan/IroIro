//----------------------------------------------------------------------
// staples.h
//----------------------------------------------------------------------
#ifndef STAPLES_INCLUDED
#define STAPLES_INCLUDED

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

class Staples {
private:
  int Nc_;
  int Ndim_;
  int Nvol_, Lvol_;
  const Format::Format_G& gf_;

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
     format1d_(new Format::Format_G(Nvol_,1)){}

  ~Staples(){delete format1d_;}

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
