/*! @file foprHermFactory.hpp
 *  @brief factory of fopr_Herm
 */
#ifndef FOPRHERMFACTORY_INCLUDED
#define FOPRHERMFACTORY_INCLUDED

#include "fopr.h"
#include "fopr_chebyshev_DdagD.h"
#include "pugi_interface.h"
#include "Dirac_ops/dirac_Operator_Factory.hpp"

class FoprHermFactory{
public:
  virtual Fopr_Herm* getFoprHerm(const DiracWilsonLike*) = 0;
  virtual ~FoprHermFactory(){}
};

class FoprHermFactory_H: public FoprHermFactory{
public:
  FoprHermFactory_H(XML::node node){}
  Fopr_Herm* getFoprHerm(const DiracWilsonLike* D){ return new Fopr_H(D); }
};

class FoprHermFactory_DdagD: public FoprHermFactory{
public:
  FoprHermFactory_DdagD(XML::node node){}
  Fopr_Herm* getFoprHerm(const DiracWilsonLike* D){ return new Fopr_DdagD(D); }
};

class FoprHermFactory_DDdag: public FoprHermFactory{
public:
  FoprHermFactory_DDdag(XML::node node){}
  Fopr_Herm* getFoprHerm(const DiracWilsonLike* D){ return new Fopr_DDdag(D); }
};

class FoprHermFactory_Chebyshev_DdagD: public FoprHermFactory{
  XML::node fopr_node_;
public:
  FoprHermFactory_Chebyshev_DdagD(XML::node node):fopr_node_(node){}
  Fopr_Herm* getFoprHerm(const DiracWilsonLike* D){ 
    return new Fopr_Chebyshev_DdagD(fopr_node_,D);
  }
};

namespace FoprHerm{ 
  FoprHermFactory* createFoprHermFactory(XML::node);
}
#endif
