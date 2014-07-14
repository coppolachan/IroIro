/*!@file chebyshevAccelFuncFactory.hpp
 * @brief factory of fopr_Chebyshev whose kernel is DdagDLin
 This is especially used for the chebyshev acceleration in the eigenmodes calc.
 */
#ifndef CHEBYSHEVACCELFUNCFACTORY_INCLUDED
#define CHEBYSHEVACCELFUNCFACTORY_INCLUDED

#include "eigenCalcGeneral.hpp"
#include "Fopr/fopr_Linear.h"
#include "Fopr/fopr_QuadLinear.h"
#include "Fopr/fopr_Chebyshev.h"

class ChebyshevAccelFuncFactory: public OpFunc{
  //// internal classes
  class LinearFuncFactory: public OpFunc{
    double slp_,icpt_;
  public:
    LinearFuncFactory(double slp,double icpt):slp_(slp),icpt_(icpt){}
    Fopr_Linear* getOp(const Fopr_Herm* op)const{ 
      return new Fopr_Linear(slp_,icpt_,op); 
    }
  };
  
  class QuadLinearFuncFactory: public OpFunc{
    double slp_,icpt_;
  public:
    QuadLinearFuncFactory(double slp,double icpt):slp_(slp),icpt_(icpt){}
    Fopr_QuadLinear* getOp(const Fopr_Herm* op)const{ 
      return new Fopr_QuadLinear(slp_,icpt_,op); 
    }
  };

  ////
  std::auto_ptr<OpFunc> mapFact_;  
  mutable std::auto_ptr<Fopr_Herm> mapOp_;  
  XML::node cbNode_;  /// node to chebyshev parameters
public:
  ChebyshevAccelFuncFactory(XML::node setupNode);
  Fopr_Chebyshev* getOp(const Fopr_Herm* op)const;
};

#endif
