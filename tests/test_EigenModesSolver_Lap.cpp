/*! @file test_EigenModesSolver_Lap.cpp
 * @brief implementation of the Test_EigenModesSolver_Lap class
 */
#include "test_EigenModesSolver_Lap.hpp"
#include "EigenModes/eigenCalcGeneral.hpp"

#include "format_S.h"

using namespace std;

int Test_EigenModesSolver_Lap::run(){
  XML::node eig_node= input_.node;
  XML::descend(eig_node,"EigenModesCalc");
  EigenCalcGeneral eigen(eig_node);

  try{
    eigen.do_calc(input_.config);
    eigen.output_bin3D<Format::Format_S>(input_.output);
  
  }
  catch(const char* error){  CCIO::cout<<error<<"\n"; }
  return 0;
}
