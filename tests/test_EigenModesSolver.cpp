/*! @file test_EigenModesSolver.cpp
 * @brief implementation of the Test_EigenModesSolver class
 */
#include "test_EigenModesSolver.hpp"
#include "EigenModes/eigenCalcGeneral.hpp"

using namespace std;

int Test_EigenModesSolver::run(){
  
  XML::node eig_node= input_.node;
  XML::descend(eig_node,"EigenModesCalc");
  EigenCalcGeneral eigen(eig_node);
  
  try{
    eigen.do_calc(input_.gconf);
    eigen.output_bin(input_.output);
  }catch(const char* error){
    CCIO::cout<<error<<"\n";
  }
  return 0;
}
