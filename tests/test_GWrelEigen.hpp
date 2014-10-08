/*! @file test_GWrelEigen.hpp
 *  @brief A study of the GW-relation via the Dirac spectrum
 */
#ifndef TEST_GWRELEIGEN_INCLUDED
#define TEST_GWRELEIGEN_INCLUDED

#include "Dirac_ops/dirac_Operator_FactoryCreator.hpp"
#include <vector>
#include "include/iroiro_code.hpp"
#include "Communicator/comm_io.hpp"
#include "tests.hpp"

class Field;

class Test_GWrelEigen{
  const Measurements::Input input_;
 public:
  Test_GWrelEigen(const Measurements::Input& input):input_(input){
    CCIO::cout<<"Test_GWrelEigen called\n";    
  }
  int run();
};

#endif
