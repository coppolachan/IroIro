/*! @file foprHermFactory.hpp
 *  @brief factory of fopr_Herm
 */
#ifndef FOPRHERMFACTORY_INCLUDED
#define FOPRHERMFACTORY_INCLUDED

#include "fopr.h"
#include "pugi_interface.h"
#include "Dirac_ops/dirac_Operator_Factory.hpp"
#include "Dirac_ops/dirac_Operator_FactoryCreator.hpp"
#include "Scalar_ops/scalar_Operator_Factory.hpp"
#include "Scalar_ops/scalar_Operator_FactoryCreator.hpp"
#include <memory>

class FoprHermFactory{
public:
  virtual Fopr_Herm* getFoprHerm(InputConfig&) = 0;
  virtual ~FoprHermFactory(){}
};

class FoprHermFactory_HD: public FoprHermFactory{
  std::auto_ptr<DiracWilsonLikeFactory> Dfact_;
  std::auto_ptr<DiracWilsonLike> D_;
public:
  FoprHermFactory_HD(XML::node node){
    XML::descend(node,"WilsonLikeDirac",MANDATORY);
    Dfact_.reset(Diracs::createDiracWilsonLikeFactory(node));
  }
  Fopr_HD* getFoprHerm(InputConfig& input){ 
    D_.reset(Dfact_->getDirac(input));
    return new Fopr_HD(D_.get()); 
  }
};

class FoprHermFactory_H: public FoprHermFactory{
  std::auto_ptr<DiracWilsonLikeFactory> Dfact_;
  std::auto_ptr<DiracWilsonLike> D_;
public:
  FoprHermFactory_H(XML::node node){
    XML::descend(node,"WilsonLikeDirac",MANDATORY);
    Dfact_.reset(Diracs::createDiracWilsonLikeFactory(node));
  }
  Fopr_H* getFoprHerm(InputConfig& input){ 
    D_.reset(Dfact_->getDirac(input));
    return new Fopr_H(D_.get()); 
  }
};

class FoprHermFactory_H5d: public FoprHermFactory{
  std::auto_ptr<DiracDomainWall5dFactory> Dfact_;
  std::auto_ptr<Dirac_DomainWall> D_;
public:
  FoprHermFactory_H5d(XML::node node){
    XML::descend(node,"WilsonLikeDirac",MANDATORY);
    Dfact_.reset(new DiracDomainWall5dFactory(node));
  }
  Fopr_H5d* getFoprHerm(InputConfig& input){ 
    D_.reset(Dfact_->getDirac(input));
    return new Fopr_H5d(D_.get()); 
  }
};

class FoprHermFactory_DdagD: public FoprHermFactory{
  std::auto_ptr<DiracWilsonLikeFactory> Dfact_;
  std::auto_ptr<DiracWilsonLike> D_;
public:
  FoprHermFactory_DdagD(XML::node node){
    XML::descend(node,"WilsonLikeDirac",MANDATORY);
    Dfact_.reset(Diracs::createDiracWilsonLikeFactory(node));
  }
  Fopr_DdagD* getFoprHerm(InputConfig& input){ 
    D_.reset(Dfact_->getDirac(input));
    return new Fopr_DdagD(D_.get()); 
  }
};

class FoprHermFactory_DDdag: public FoprHermFactory{
  std::auto_ptr<DiracWilsonLikeFactory> Dfact_;
  std::auto_ptr<DiracWilsonLike> D_;
public:
  FoprHermFactory_DDdag(XML::node node){
    XML::descend(node,"WilsonLikeDirac",MANDATORY);
    Dfact_.reset(Diracs::createDiracWilsonLikeFactory(node));
  }
  Fopr_DDdag* getFoprHerm(InputConfig& input){ 
    D_.reset(Dfact_->getDirac(input));
    return new Fopr_DDdag(D_.get()); 
  }
};

class FoprHermFactory_Scalar: public FoprHermFactory{
  std::auto_ptr<ScalarOpFactory> Sfact_;
  std::auto_ptr<ScalarOp> S_;
public:
  FoprHermFactory_Scalar(XML::node node){
    XML::descend(node,"ScalarOp",MANDATORY);
    assert(ScOps::createScalarOpFactory(node));
    Sfact_.reset(ScOps::createScalarOpFactory(node));
  }
  Fopr_Scalar* getFoprHerm(InputConfig& input){ 
    S_.reset(Sfact_->getSop(input));
    return new Fopr_Scalar(S_.get()); 
  }
};

#endif
