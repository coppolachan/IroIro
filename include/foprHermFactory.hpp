/*! @file foprHermFactory.hpp
 *  @brief factory of fopr_Herm
 */
#ifndef FOPRHERMFACTORY_INCLUDED
#define FOPRHERMFACTORY_INCLUDED

#include "fopr.h"
#include "fopr_Chebyshev.h"
#include "fopr_DdagDLinear.h"
#include "pugi_interface.h"
#include "Dirac_ops/dirac_Operator_Factory.hpp"

class FoprHermFactory{
public:
  virtual Fopr_Herm* getFoprHerm(const DiracWilsonLike*) = 0;
  virtual ~FoprHermFactory(){}
};

class FoprHermFactory_H: public FoprHermFactory{
public:
  Fopr_H* getFoprHerm(const DiracWilsonLike* D){ return new Fopr_H(D); }
};

class FoprHermFactory_DdagD: public FoprHermFactory{
public:
  Fopr_DdagD* getFoprHerm(const DiracWilsonLike* D){ return new Fopr_DdagD(D);}
};

class FoprHermFactory_DDdag: public FoprHermFactory{
public:
  Fopr_DDdag* getFoprHerm(const DiracWilsonLike* D){ return new Fopr_DDdag(D);}
};

class FoprHermFactory_DdagDLinear: public FoprHermFactory{
  double slp_,icp_;
public:
  FoprHermFactory_DdagDLinear():slp_(1.0),icp_(0.0){}
  FoprHermFactory_DdagDLinear(double slope,double intercept)
    :slp_(slope),icp_(intercept){}
  
  Fopr_DdagDLinear* getFoprHerm(const DiracWilsonLike* D){ 
    return new Fopr_DdagDLinear(D,slp_,icp_); 
  }
};

class FoprHermFactory_Chebyshev: public FoprHermFactory{
  const Fopr_Herm* kernel_;
  int N_;
  FoprHermFactory_Chebyshev(FoprHermFactory_Chebyshev&);
  FoprHermFactory_Chebyshev& operator=(FoprHermFactory_Chebyshev&);
public:
  FoprHermFactory_Chebyshev(const Fopr_Herm* kernel,int Npoly)
    :kernel_(kernel),N_(Npoly){}
  
  Fopr_Chebyshev* getFoprHerm(const DiracWilsonLike* D){ 
    return new Fopr_Chebyshev(kernel_,N_);
  }
};

#endif
