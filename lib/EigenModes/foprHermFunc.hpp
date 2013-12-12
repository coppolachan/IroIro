#ifndef FOPRHERMFUNC_INCLUDED
#define FOPRHERMFUNC_INCLUDED

#include "include/fopr_Chebyshev.h"
#include "include/fopr_Linear.h"
#include "include/fopr_QuadLinear.h"
#include <iostream>
#include <memory>

class FoprHermFunc{
public:
  virtual Fopr_Herm* getFoprHerm(const Fopr_Herm*)const = 0;
  virtual ~FoprHermFunc(){}
};

class FoprNULLfunc: public FoprHermFunc{
public:
  Fopr_Herm* getFoprHerm(const Fopr_Herm* foprH)const { return NULL;}
};

class FoprChebyshevFunc: public FoprHermFunc{
  XML::node cbnode_;
public:
  FoprChebyshevFunc(XML::node node):cbnode_(node){}

  Fopr_Chebyshev* getFoprHerm(const Fopr_Herm* foprH)const {
    return new Fopr_Chebyshev(cbnode_,foprH);
  }
};

class FoprLinearFunc: public FoprHermFunc{
  XML::node lnode_;
  double slp_;
  double icpt_;
public:
  FoprLinearFunc(XML::node node):lnode_(node){}
  FoprLinearFunc(double slp,double icpt)
    :slp_(slp),icpt_(icpt),lnode_(NULL){}
  
  Fopr_Linear* getFoprHerm(const Fopr_Herm* foprH)const{ 
    if(lnode_!=NULL) return new Fopr_Linear(lnode_,foprH); 
    else             return new Fopr_Linear(slp_,icpt_,foprH); 
  }
};

class FoprQuadLinearFunc: public FoprHermFunc{
  XML::node lnode_;
  double slp_;
  double icpt_;
public:
  FoprQuadLinearFunc(XML::node node):lnode_(node){}
  FoprQuadLinearFunc(double slp,double icpt)
    :slp_(slp),icpt_(icpt),lnode_(NULL){}
  
  Fopr_QuadLinear* getFoprHerm(const Fopr_Herm* foprH)const{ 
    if(lnode_!=NULL) return new Fopr_QuadLinear(lnode_,foprH); 
    else             return new Fopr_QuadLinear(slp_,icpt_,foprH); 
  }
};

#endif
