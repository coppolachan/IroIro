//------------------------------------------------------------------------
// test_EigenModes_IRL.h
//------------------------------------------------------------------------
#ifndef TEST_EIGENSOLVER_INCLUDED
#define TEST_EIGENSOLVER_INCLUDED

#include "include/common_code.hpp"
#include "tests/tests.hpp"

class SortEigen_low;
class SortEigen_high;

class Test_EigenModes_IRL: public TestGeneral{
 private:
  XML::node Eigen_node;
  GaugeField& u_;
 public:
  Test_EigenModes_IRL(XML::node node,
		      GaugeField& conf):Eigen_node(node),
					u_(conf){}
  int lowlying();
  int highest();
  int chebyshev();

  int run();
};

#endif
