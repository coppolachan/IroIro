#include "field.h"
#include "common_fields.hpp"

class Laplacian {
  GaugeField *u;
  
public:
  Laplacian(GaugeField*);
  
  FermionField1sp apply(const FermionField1sp&);

};
