//------------------------------------------------------------------------
// test_EigenModes_IRL.h
//------------------------------------------------------------------------
#ifndef TEST_EIGENSOLVER_INCLUDED
#define TEST_EIGENSOLVER_INCLUDED

#include "include/common_code.hpp"


class SortEigen_low;
class SortEigen_high;

class Test_EigenModes_IRL{
 private:
  Field& u_;
 public:
  Test_EigenModes_IRL(Field& conf):u_(conf){}
  int lowlying();
  int highest();
  int chebyshev();

  int run();
};

#endif
