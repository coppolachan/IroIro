/*!@file quark_prop_meas_factory.cpp
 * @brief implementations quark_propagator factories
 */
#include "quark_prop_meas_factory.hpp"
#include "include/common_fields.hpp"
#include "Dirac_ops/dirac_Operator_FactoryCreator.hpp"
#include "PugiXML/xmlUtilities.hpp"
#include <string.h>
#include <sstream>

using namespace std;

/// QPropFactory,   some leaking expected - check!
QPropFactory::QPropFactory(const XML::node& node):node_(node){ 
  XML::descend(node_,"Kernel");
  Dfactory_.save(Diracs::createDiracWilsonLikeFactory(node_));
  XML::next_sibling(node_,"Solver");
  slvFactory_.save(Solvers::createSolverFactory(node_));
}

QPropFactory::QPropFactory(const XML::node& node,double mass):node_(node){ 
  XML::descend(node_,"Kernel");

  stringstream mass_val;
  mass_val << mass;
  if(! (node_.child("mass").first_child().set_value(mass_val.str().c_str()))){ 
    CCIO::cout<<"QPropFactory: set_value failed\n";
    abort();
  }
  Dfactory_.save(Diracs::createDiracWilsonLikeFactory(node_));
  XML::next_sibling(node_,"Solver");
  slvFactory_.save(Solvers::createSolverFactory(node_));
}

QuarkPropagator* QPropFactory::createQuarkProp(const InputConfig& input){
  Kernel_.save(Dfactory_.get()->getDirac(input));
  DdagD_.save(new Fopr_DdagD(Kernel_.get()));
  Solv_.save(slvFactory_.get()->getSolver(DdagD_.get()));
  return new Qprop(Kernel_.get(),Solv_.get());
}

/// QProp_Factory_EvenOdd 
QPropFactory_EvenOdd::QPropFactory_EvenOdd(const XML::node& node):node_(node){
  XML::descend(node_,"Kernel");
  Dfactory_.save(new DiracWilsonEvenOddFactory(node_));
  XML::next_sibling(node_,"Solver");
  slvFactory_.save(Solvers::createSolverFactory(node_));
}

QPropFactory_EvenOdd::QPropFactory_EvenOdd(const XML::node& node,double mass)
  :node_(node){
  XML::descend(node_,"Kernel");

  stringstream mass_val;
  mass_val << mass;
  if(! (node_.child("mass").first_child().set_value(mass_val.str().c_str()))){ 
    CCIO::cout<<"QPropFactory_EvenOdd: set_value failed\n";
    abort();
  }
  Dfactory_.save(new DiracWilsonEvenOddFactory(node_));
  XML::next_sibling(node_,"Solver");
  slvFactory_.save(Solvers::createSolverFactory(node_));
}

QuarkPropagator* QPropFactory_EvenOdd::createQuarkProp(const InputConfig& input){
  Kernel_.save(Dfactory_.get()->getDirac(input));
  DdagD_.save(new Fopr_DdagD(Kernel_.get()));
  Solv_.save(slvFactory_.get()->getSolver(DdagD_.get()));
  return new Qprop_EvenOdd(Kernel_.get(),Solv_.get());
}

/// QPropDWFFactory
QPropDWFFactory::QPropDWFFactory(const XML::node& node):node_(node){
  XML::descend(node_,"Kernel4d");
  Dfactory_.save(Diracs::createDiracDWF4dFactory(node_));
}

QPropDWFFactory::QPropDWFFactory(const XML::node& node,double mass):node_(node){
  XML::descend(node_,"Kernel4d");

  stringstream mass_val;
  mass_val << mass;
  if(! (node_.child("Kernel5d").child("mass").first_child().set_value(mass_val.str().c_str()))){ 
    CCIO::cout<<"QPropDWFFactory: set_value failed\n";
    abort();
  }
  Dfactory_.save(Diracs::createDiracDWF4dFactory(node_));

  double mq;
  XML::descend(node_,"Kernel5d");
  XML::read(node_,"mass",mq,MANDATORY);
  CCIO::cout<<"QPropDWFFactory: mass= "<<mq<<"\n";
}

QpropDWF* QPropDWFFactory::createQuarkProp(const InputConfig& input){
  DWF4D_.save(Dfactory_.get()->getDirac(input));
  return new QpropDWF(*DWF4D_.get());
}

/////////////////// factory creator  ///////////////////////

namespace QuarkPropagators {
  QuarkPropagatorFactory* createQuarkPropagatorFactory(const XML::node& node){
    
    XML::nullCheck(node,"QuarkPropagator");
    const char* qpname = node.attribute("name").value();
    
    if(!strcmp(qpname,"Qprop"))        return new QPropFactory(node);
    if(!strcmp(qpname,"Qprop_EvenOdd"))return new QPropFactory_EvenOdd(node);
    if(!strcmp(qpname,"QpropDWF"))     return new QPropDWFFactory(node);
    
    XML::stopMsg(node,"QuarkPropagator");
  }  
  
  QuarkPropagatorFactory* createQuarkPropagatorFactory(const XML::node& node,
						       double mass){
    XML::nullCheck(node,"QuarkPropagator");
    const char* qpname = node.attribute("name").value();
    
    if(!strcmp(qpname,"Qprop"))        return new QPropFactory(node,mass);
    if(!strcmp(qpname,"Qprop_EvenOdd"))return new QPropFactory_EvenOdd(node,mass);
    if(!strcmp(qpname,"QpropDWF"))     return new QPropDWFFactory(node,mass);

    XML::stopMsg(node,"QuarkPropagator");
  }  
}
