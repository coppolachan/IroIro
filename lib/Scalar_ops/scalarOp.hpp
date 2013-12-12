#ifndef SCALAROP_INCLUDED
#define SCALAROP_INCLUDED

#include "include/field.h"

class ScalarOp{
public:
  virtual ~ScalarOp(){}
  virtual size_t fsize() const =0;
  virtual const Field mult(const Field&)const = 0;
};
#endif
