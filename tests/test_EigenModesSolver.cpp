/*! @file test_EigenModesSolver.cpp
 * @brief implementation of the Test_EigenModesSolver class
 */
#include "test_EigenModesSolver.hpp"
#include "EigenModes/eigenCalcGeneral.hpp"

#include "Dirac_ops/BoundaryConditions/boundaryCond.hpp"

#include "format_F.h"

using namespace std;

int Test_EigenModesSolver::run(){
  BoundaryCond* BC;  
  XML::node eig_node= input_.node;
  XML::descend(eig_node,"EigenModesCalc");
  EigenCalcGeneral eigen(eig_node);

  //Apply boundary condition                                                           
  bool AntiPeriodic = false; // default                                             
  XML::read(eig_node, "AntiPeriodicBC", AntiPeriodic);
  
  try{
    if (AntiPeriodic){
      CCIO::cout << "Appying Antiperiodic BC\n";
      BC = new BoundaryCond_antiPeriodic(TDIR);
      BC->apply_bc(*input_.gconf);
    }
    
    InputConfig input = input_.getConfig();
    eigen.do_calc(input);
    eigen.output_bin<Format::Format_F>(input_.output);
  
  }
  catch(const char* error){
    CCIO::cout<<error<<"\n";
  }

  return 0;
}
