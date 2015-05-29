/*!
 * @file test_Solver_HDCG.hpp
 *
 * @brief Declaration of classes for testing the BFM HDCG classes
 *
 */
#ifndef TEST_SOLVER_BFM_HDCG_INCLUDED
#define TEST_SOLVER_BFM_HDCG_INCLUDED

#include "include/iroiro_code.hpp"
#include "tests/tests.hpp"
#include "Dirac_ops/BFM_Wrapper/dirac_BFM_HDCG.hpp"

class Test_Solver_HDCG: public TestGeneral{
private:
  const char* test_name;
  XML::node DWFnode;
  XML::node top_node;
  GaugeField conf_;


public:
  Test_Solver_HDCG(XML::node node,GaugeField conf):DWFnode(node),
						  conf_(conf){
    test_name = "Test_HDCG";
    top_node = DWFnode; 
    XML::descend(DWFnode, test_name, MANDATORY);
  }

  int run();
  int launch_test(FermionField&, FermionField&, Dirac_BFM_HDCG_Wrapper*);

};

#endif
