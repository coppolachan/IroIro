/*! @file chiral_condensate_abs.hpp
 * @brief Declaration of Chiralc condensate virtual base class */
#ifndef CHIRCOND_H_
#define CHIRCOND_H_

#include "source.hpp"

/*! @class ChiralCondensate
  * @brief Calculates the Chiral condensate  \f$ \rangle \bar\psi \psi \f$ */

class ChiralCondensate {
public:
  /*! @brief Public destructor */
  virtual ~ChiralCondensate(){}
  
  /*! @brief Calculates the condensate */
  virtual double calc(Source& src) const = 0;
};
#endif

