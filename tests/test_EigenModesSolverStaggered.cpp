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
    eigen.do_calc(input_.config);
    eigen.output_bin<Format::Format_S>(input_.output,input_.app_out);
  }
  catch(const char* error){
    CCIO::cout<<error<<"\n";
  }

  return 0;
}
