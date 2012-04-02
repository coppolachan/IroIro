/*! 
  @file findminmax.hpp

  @brief Find the minimum and maximum eigenmode of a real matrix

*/

#include <memory>

#include "include/fopr.h"
#include "Tools/randNum.h"


class findMinMax {
  Fopr* Kernel_;
  std::auto_ptr<const RandNum> rand_;
  size_t vect_size;

  findMinMax(); //hide default constructor
public:
  /*! @brief Default constructor - Fopr should have real eigenmodes */
  explicit findMinMax(Fopr*, const RandNum*, size_t) ;

  ~findMinMax();
  double findMin() const;
  double findMax() const;
};
