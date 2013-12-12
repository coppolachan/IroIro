#ifndef SCALAR_OPERATOR_FACTORY_INCLUDED
#define SCALAR_OPERATOR_FACTORY_INCLUDED

#include "scalarOp.hpp"
#include "laplacian.hpp"
#include "pugi_interface.h"
#include "include/inputConfig.hpp"
#include <string.h>

class ScalarOpFactory{
public:
  virtual ScalarOp* getSop(InputConfig& )const = 0;
  virtual ~ScalarOpFactory(){}
};

class LaplacianFactory: public ScalarOpFactory{
  XML::node lnode_;
public:
  LaplacianFactory(const XML::node& node):lnode_(node){}
  Laplacian* getSop(InputConfig& input)const{
    return new Laplacian(lnode_,input.getGconf());
  }
};

#endif
