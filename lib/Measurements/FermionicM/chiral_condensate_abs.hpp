/*! @file chiral_condensate_abs.hpp
 * @brief Declaration of Chiralc condensate virtual base class */
#ifndef CHIRCOND_H_
#define CHIRCOND_H_

#include "source.hpp"
#include "include/field.h"

/*!
 * @class ChiralCondensate
 * @brief Calculates the Chiral condensate  \f$ \rangle \bar\psi \psi \f$ 
*/

class ChiralCondensate {
  virtual Field invert (Field&) const = 0;
  virtual Field gamma5 (Field&) const = 0;
  virtual int fsize()const = 0;
public:
  /*! @brief Public destructor */
  virtual ~ChiralCondensate(){}
  
  /*! @brief Calculates the condensate */
  double calc(Source& src, const int) const;
  double connected_susc(Source& src, const int) const;
};
#endif

