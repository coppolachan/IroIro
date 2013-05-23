/*!  @file quark_prop_meas_factory.hpp 
  @brief QuarkPropagators measurements operators factories
*/
#ifndef QUARKPROP_MEAS_FACT_
#define QUARKPROP_MEAS_FACT_

#include "fermion_meas_factory_abs.hpp"

#include "include/pugi_interface.h"
#include "Tools/RAIIFactory.hpp"

// Concrete objects definitions
#include "qprop.hpp"
#include "qprop_EvenOdd.hpp"
#include "qprop_DomainWall.hpp"
#include "Dirac_ops/dirac_Operator_FactoryCreator.hpp"
#include "Solver/solver_Factory.hpp"

/////////////////////////////////////////////////
/*! @brief Concrete class for creating Quark Propagator Qprop operator */

class QPropFactory : public QuarkPropagatorFactory {
  RaiiFactoryObj<DiracWilsonLikeFactory> DiracObj;
  RaiiFactoryObj<SolverFactory> SolverObj;

  RaiiFactoryObj<DiracWilsonLike> Kernel;
  RaiiFactoryObj<Solver> Solv;
public:
  QPropFactory(XML::node node){
    //some leaking expected - check!
    XML::descend(node,"Kernel");
    DiracObj.save(Diracs::createDiracWilsonLikeFactory(node));
    XML::next_sibling(node,"Solver");
    SolverObj.save(Solvers::createSolverFactory(node));
  }

  QuarkPropagator* getQuarkProp(GaugeField& Field){
    Kernel.save(DiracObj.get()->getDiracWL(&(Field.data)));
    Solv.save(SolverObj.get()->getSolver(new Fopr_DdagD(Kernel.get())));

    return new Qprop(Kernel.get(), Solv.get());
  }
};

/*!
 @brief Concrete class for creating Quark Propagator Qprop_EvenOdd operator
*/
class QPropFactory_EvenOdd : public QuarkPropagatorFactory {
  RaiiFactoryObj<DiracWilsonLikeEvenOddFactory> DiracObj;
  RaiiFactoryObj<SolverFactory> SolverObj;

  RaiiFactoryObj<DiracWilsonLike_EvenOdd> Kernel;
  RaiiFactoryObj<Solver> Solv;
public:
  QPropFactory_EvenOdd(XML::node node){
    XML::descend(node,"Kernel");
    DiracObj.save(new DiracWilsonEvenOddFactory(node));
    XML::next_sibling(node,"Solver");
    SolverObj.save(Solvers::createSolverFactory(node));
  }

  QuarkPropagator* getQuarkProp(GaugeField& Field){
    Kernel.save(DiracObj.get()->getDiracEO(&(Field.data)));
    Solv.save(SolverObj.get()->getSolver(new Fopr_DdagD(Kernel.get())));
    return new Qprop_EvenOdd(Kernel.get(), Solv.get());
  }
};

/*!
 @brief Concrete class for creating Quark Propagator QpropDWF operator
*/
class QPropDWFFactory : public QuarkPropagatorFactory {
  RaiiFactoryObj<DiracDWF4dFactory> DiracFactory_;
  RaiiFactoryObj<Dirac_optimalDomainWall_4D> DWF4D_;

public:
  QPropDWFFactory(XML::node node){
    XML::descend(node,"Kernel4d");
    DiracFactory_.save(Diracs::createDiracDWF4dFactory(node));
  }

  QuarkPropagator* getQuarkProp(GaugeField& Field){
    return getQuarkPropDW(Field);
  }

  QpropDWF* getQuarkPropDW(GaugeField& Field){
    DWF4D_.save(DiracFactory_.get()->getDirac4D(&(Field.data)));
    return new QpropDWF(*DWF4D_.get());
  }
};

#endif 
