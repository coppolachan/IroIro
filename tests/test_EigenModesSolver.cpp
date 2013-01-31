/*! @file test_EigenModesSolver.cpp
 * @brief implementation of the Test_EigenModesSolver class
 */
#include "test_EigenModesSolver.hpp"
#include "EigenModes/eigenCalcGeneral.hpp"

using namespace std;

int Test_EigenModesSolver::run(){
  XML::node eig_node= node_;
  XML::descend(eig_node,"EigenModesCalc");
  EigenCalcGeneral eigen(eig_node);
  
  try{
    eigen.do_calc(&(conf_.data));
    eigen.output_bin(output_);
  }catch(const char* error){
    CCIO::cout<<error<<"\n";
  }
  return 0;
}
