/*! @file test_EigenModesSolver.cpp
 * @brief implementation of the Test_EigenModesSolver class
 */
#include "test_EigenModesSolver.hpp"
#include "EigenModes/eigenCalcGeneral.hpp"
#include "format_F.h"

using namespace std;

int Test_EigenModesSolver::run(){
  
  XML::node eig_node= input_.node;
  XML::descend(eig_node,"EigenModesCalc");
  EigenCalcGeneral eigen(eig_node);
  
  try{
    InputConfig input = input_.getConfig();
    eigen.do_calc(input);
    //eigen.output_bin<Format::Format_F>(input_.output);
    eigen.output_bin3D<Format::Format_S>(input_.output);
  }catch(const char* error){
    CCIO::cout<<error<<"\n";
  }
  return 0;
}
