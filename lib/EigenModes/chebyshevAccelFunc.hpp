/*!@file chebyshevAccelFunc.hpp
 * @brief factory of fopr_Chebyshev whose kernel is DdagDLin
 This is especially used for the chebyshev acceleration in the eigenmodes calc.
 */
#ifndef CHEBYSHEVACCELFUNC_INCLUDED
#define CHEBYSHEVACCELFUNC_INCLUDED

#include "foprHermFunc.hpp"

class ChebyshevAccelFunc: public FoprHermFunc{
  std::auto_ptr<FoprHermFunc> mapFact_;  
  mutable std::auto_ptr<Fopr_Herm> mapOp_;  

  std::auto_ptr<FoprChebyshevFunc> cbfunc_;  
public:
  ChebyshevAccelFunc(XML::node node);
  Fopr_Chebyshev* getFoprHerm(const Fopr_Herm* foprH)const;
};

#endif
