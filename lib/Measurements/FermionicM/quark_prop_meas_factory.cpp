/*!@file quark_prop_meas_factory.cpp
 * @brief implementations quark_propagator factories
 */
#include "quark_prop_meas_factory.hpp"
#include "include/common_fields.hpp"
#include "Dirac_ops/dirac_Operator_FactoryCreator.hpp"
#include <string.h>

/// QPropFactory,   some leaking expected - check!
QPropFactory::QPropFactory(XML::node node){ 
  XML::descend(node,"Kernel");
  Dfactory_.save(Diracs::createDiracWilsonLikeFactory(node));
  XML::next_sibling(node,"Solver");
  slvFactory_.save(Solvers::createSolverFactory(node));
}

QuarkPropagator* QPropFactory::createQuarkProp(InputConfig& input){
  Kernel_.save(Dfactory_.get()->getDirac(input));
  Solv_.save(slvFactory_.get()->getSolver(new Fopr_DdagD(Kernel_.get())));
  return new Qprop(Kernel_.get(),Solv_.get());
}

/// QProp_Factory_EvenOdd 
QPropFactory_EvenOdd::QPropFactory_EvenOdd(XML::node node){
  XML::descend(node,"Kernel");
  Dfactory_.save(new DiracWilsonEvenOddFactory(node));
  XML::next_sibling(node,"Solver");
  slvFactory_.save(Solvers::createSolverFactory(node));
}

QuarkPropagator* QPropFactory_EvenOdd::createQuarkProp(InputConfig& input){
  Kernel_.save(Dfactory_.get()->getDirac(input));
  Solv_.save(slvFactory_.get()->getSolver(new Fopr_DdagD(Kernel_.get())));
  return new Qprop_EvenOdd(Kernel_.get(),Solv_.get());
}

/// QPropDWFFactory
QPropDWFFactory::QPropDWFFactory(XML::node node){
  XML::descend(node,"Kernel4d");
  Dfactory_.save(Diracs::createDiracDWF4dFactory(node));
}

QpropDWF* QPropDWFFactory::createQuarkProp(InputConfig& input){
  DWF4D_.save(Dfactory_.get()->getDirac(input));
  return new QpropDWF(*DWF4D_.get());
}

/////////////////// factory creator  ///////////////////////

namespace QuarkPropagators {
  QuarkPropagatorFactory* createQuarkPropagatorFactory(XML::node node) {
    
    if(node !=NULL){
      const char* qprop_name = node.attribute("name").value();
        
      if(!strcmp(qprop_name, "Qprop"))  
        return new QPropFactory(node);
      if(!strcmp(qprop_name, "Qprop_EvenOdd")) 
        return new QPropFactory_EvenOdd(node);
      if(!strcmp(qprop_name, "QpropDWF")) 
        return new QPropDWFFactory(node);

      CCIO::cerr << "No Quark Propagator available with name ["
		 << qprop_name << "]" << std::endl;
      abort();
    }else{
      CCIO::cout << "Requested node is missing in input file "
		 << "(QuarkPropagator Object)\n" 
		 << "Request by " << node.parent().name() << std::endl;
      abort();
    }
  }  
}
