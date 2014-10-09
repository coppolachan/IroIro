#ifndef SCALAR_OPERATOR_FACTORY_INCLUDED
#define SCALAR_OPERATOR_FACTORY_INCLUDED

#include "scalarOp.hpp"
#include "laplacian.hpp"
#include "laplacian4Ds.hpp"
#include "pugi_interface.h"
#include "include/inputConfig.hpp"
#include <string.h>

class ScalarOpFactory{
public:
  virtual ScalarOp* getSop(const InputConfig& )const = 0;
  virtual ~ScalarOpFactory(){}
};

class LaplacianFactory: public ScalarOpFactory{
  XML::node lnode_;
public:
  LaplacianFactory(const XML::node& node):lnode_(node){}
  Laplacian* getSop(const InputConfig& input)const{
    return new Laplacian(lnode_,input.getGconf());
  }
};

class Laplacian4DFfactory: public ScalarOpFactory{
public:
  Laplacian4DFfactory(){}
  Laplacian4DF* getSop(const InputConfig& input)const{
    CCIO::cout<<"Laplacian4DFfactory called\n";
    return new Laplacian4DF(input.getGconf());
  }
};

class Laplacian4DSfactory: public ScalarOpFactory{
public:
  Laplacian4DSfactory(){}
  Laplacian4DS* getSop(const InputConfig& input)const{
    return new Laplacian4DS(input.getGconf());
  }
};

#endif
