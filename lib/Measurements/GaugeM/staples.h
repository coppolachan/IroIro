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
  Format::Format_G* sf_;

  // for variables with the direction unfixed 
  SUNmat u(const Field& g,int site,int dir) const;
  SUNmat u_dag(const Field& g,int site,int dir) const;

  // for variables with a specific direction
  SUNmat u(const Field& g,int site) const;
  SUNmat u(const std::valarray<double>& vu,int site) const;
  SUNmat u(const ShiftField& su,int site) const;

  SUNmat u_dag(const Field& g,int site) const;
  SUNmat u_dag(const std::valarray<double>& vu,int site) const;
  SUNmat u_dag(const ShiftField& su,int site) const;

public:
  Staples(const Format::Format_G& gf)
    :gf_(gf),
     Nc_(CommonPrms::instance()->Nc()),
     Ndim_(CommonPrms::instance()->Ndim()),
     Nvol_(CommonPrms::instance()->Nvol()),
     Lvol_(CommonPrms::instance()->Lvol()),
     com_(Communicator::instance()),
     idx_(SiteIndex::instance()),
     sf_(new Format::Format_G(Nvol_,1)){}

  ~Staples(){delete sf_;}

  // functions for Field class
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
