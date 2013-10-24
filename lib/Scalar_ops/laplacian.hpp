#ifndef LAPLACIAN_INCLUDED
#define LAPLACIAN_INCLUDED

#include "scalarOp.hpp"
#include "common_fields.hpp"
#include "Geometry/shiftField.hpp"
#include <cassert>

class Laplacian:public ScalarOp {
  const Field * const u_;
  int t_;
  int tsl_size_;
  Format::Format_G gfmt_;
  Format::Format_S sfmt_;

  SiteMap::Map3d map_;
  Mapping::ShiftField shiftField3d_;

  void setup();
public:
  Laplacian(int time_slice,Field* u)
    :t_(time_slice),u_(u),
     gfmt_(CommonPrms::instance()->Nvol()),
     sfmt_(CommonPrms::instance()->Nvol()/CommonPrms::instance()->Nt()){
    setup();
  }

  Laplacian(const XML::node& node,Field* u)
    :u_(u),
     gfmt_(CommonPrms::instance()->Nvol()),
     sfmt_(CommonPrms::instance()->Nvol()/CommonPrms::instance()->Nt()){
    XML::read(node,"time_slice",t_,MANDATORY);    
    setup();
  }
  
  const Field mult(const Field&)const;
  double func(double x)const{ return x;}
  size_t fsize()const {return sfmt_.size();}
};

#endif
