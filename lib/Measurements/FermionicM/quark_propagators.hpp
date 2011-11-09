/*!
 *
 * @file quark_propagators.hpp
 *
 * @brief Declaration of Quark propagator virtual base class
 *
 */

#ifndef QUARK_PROP_H_
#define QUARK_PROP_H_

#include <vector>
#include "source.h"
#include "include/field.h"

typedef std::vector<Field> prop_t;


/*!
 * @class QuarkPropagator
 *
 * @brief Calculates the Quark propagator \f$D^{-1}\f$
 *
 */

class QuarkPropagator {
public:
  /*! @brief Public destructor */
  virtual ~QuarkPropagator(){};
  
  /*! @brief Calculates the propagator */
  virtual void calc(prop_t& xq,Source& src) const = 0;
  virtual void calc(prop_t& xq,const prop_t& prp) const = 0;

};



#endif

