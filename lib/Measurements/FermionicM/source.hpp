/*!
 *
 * @file source.hpp
 *
 * @brief Definition of Source abstract base class
 *
 */

#ifndef SOURCE_INCLUDED
#define SOURCE_INCLUDED

#include "include/field.h"

/*!
 * @brief Virtual base class for sources
 *
 */
class Source{
public:
  virtual ~Source(){}
  virtual const Field mksrc(int s, int c) = 0;
  virtual const Field mksrc(const std::vector<int>& lv, int s, int c) = 0;
};



#endif
