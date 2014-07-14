/*! @file foprHermFuncFactory.hpp
 * @brief factories of functional of fopr_Herm family
 */
#ifndef FOPRHERMFUNCFACTORY_INCLUDED
#define FOPRHERMFUNCFACTORY_INCLUDED

#include "foprHermFactoryCreator.hpp"
#include "fopr_Chebyshev.h"
#include "fopr_Linear.h"
#include "fopr_QuadLinear.h"
#include "fopr_Exp.h"
#include <iostream>
#include <memory>

//// Chebyshev Function
class FoprNULLfuncFactory: public  FoprHermFactory{
public:
  Fopr_Herm* getFoprHerm(InputConfig&){return NULL; }
};

class FoprChebyshevFuncFactory: public FoprHermFactory{
  XML::node cbnode_;
  std::auto_ptr<FoprHermFactory> fhf_;
  std::auto_ptr<Fopr_Herm> hf_;
public:
  FoprChebyshevFuncFactory(XML::node node):cbnode_(node){
    XML::descend(node,"HermitianOperator", MANDATORY);
    fhf_.reset(HermiteOp::createFoprHermFactory(node));
  }

  Fopr_Chebyshev* getFoprHerm(InputConfig& input){
    hf_.reset(fhf_->getFoprHerm(input));
    return new Fopr_Chebyshev(cbnode_,hf_.get());
  }
};

//// Linear Function 
class FoprLinearFuncFactory: public FoprHermFactory{
  XML::node lnode_;
  std::auto_ptr<FoprHermFactory> fhf_;
  std::auto_ptr<Fopr_Herm> hf_;
public:
  FoprLinearFuncFactory(XML::node node):lnode_(node){
    XML::descend(node,"HermitianOperator", MANDATORY);
    fhf_.reset(HermiteOp::createFoprHermFactory(node));
  }
  Fopr_Linear* getFoprHerm(InputConfig& input){
    hf_.reset(fhf_->getFoprHerm(input));
    return new Fopr_Linear(lnode_,hf_.get()); 
  }
};

//// QuadLinear Function 
class FoprQuadLinearFuncFactory: public FoprHermFactory{
  XML::node lnode_;
  std::auto_ptr<FoprHermFactory> fhf_;
  std::auto_ptr<Fopr_Herm> hf_;
public:
  FoprQuadLinearFuncFactory(XML::node node):lnode_(node){
    XML::descend(node,"HermitianOperator", MANDATORY);
    fhf_.reset(HermiteOp::createFoprHermFactory(node));
  }
  Fopr_QuadLinear* getFoprHerm(InputConfig& input){
    hf_.reset(fhf_->getFoprHerm(input));
    return new Fopr_QuadLinear(lnode_,hf_.get()); 
  }
};

//// Exponential Function 
class FoprExpFuncFactory: public FoprHermFactory{
  XML::node enode_;
  std::auto_ptr<FoprHermFactory> fhf_;
  std::auto_ptr<Fopr_Herm> hf_;
public:
  FoprExpFuncFactory(XML::node node):enode_(node){
    XML::descend(node,"HermitianOperator", MANDATORY);
    fhf_.reset(HermiteOp::createFoprHermFactory(node));
  }
  Fopr_Exp* getFoprHerm(InputConfig& input){
    hf_.reset(fhf_->getFoprHerm(input));
    return new Fopr_Exp(enode_,hf_.get()); 
  }
};

//// LimitExponential Function 
class FoprLexpFuncFactory: public FoprHermFactory{
  XML::node enode_;
  std::auto_ptr<FoprHermFactory> fhf_;
  std::auto_ptr<Fopr_Herm> hf_;
public:
  FoprLexpFuncFactory(XML::node node):enode_(node){
    XML::descend(node,"HermitianOperator", MANDATORY);
    fhf_.reset(HermiteOp::createFoprHermFactory(node));
  }
  Fopr_Lexp* getFoprHerm(InputConfig& input){
    hf_.reset(fhf_->getFoprHerm(input));
    return new Fopr_Lexp(enode_,hf_.get()); 
  }
};

#endif
