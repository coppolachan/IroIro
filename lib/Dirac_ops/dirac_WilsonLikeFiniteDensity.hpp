/*! @file dirac_WilsonLikeFiniteDensity.hpp
 *  @brief abstract class for Wilson-Dirac operators with finite density
 */ 
#ifndef DIRACWILSONLIKEFINITEDENSITY_INCLUDED
#define DIRACWILSONLIKEFINITEDENSITY_INCLUDED

#include "dirac_WilsonLike.hpp"

class DiracWilsonLikeFiniteDensity: public DiracWilsonLike{
public:
  virtual const Field mult_Ds(const Field&)const = 0;
  virtual const Field mult_Dtp(const Field&)const = 0;
  virtual const Field mult_Dtm(const Field&)const = 0;

  virtual const Field mult_Ex(const Field& f)const{
    Field w(f.size());  // null field
    return w;
  }

  virtual ~DiracWilsonLikeFiniteDensity(){}
};
#endif
