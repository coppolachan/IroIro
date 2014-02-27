#include "test_EigModesNum.hpp"
#include "EigenModes/eigModesNum.hpp"
#include "Dirac_ops/dirac_Operator_FactoryCreator.hpp"
#include <memory>

int Test_EigModesNum::run(){
  CCIO::header(" ---- Calculating Eigen Modes Number\n");
  
  XML::node enode = input_.node;
  XML::descend(enode,"Dirac",MANDATORY);
  std::auto_ptr<DiracFactory> dfact(Diracs::createGeneralDiracFactory(enode));
  
  XML::next_sibling(enode,"EigModesNum",MANDATORY);
  InputConfig iconf = input_.getConfig();
  std::auto_ptr<Dirac> dirac(dfact.get()->getDirac(iconf));
  EigModesNum nu(enode,dirac.get(),input_.rng);
  nu.do_count();
}
