/*!
 * @file fermion_meas_factory.hpp 
 * @brief Fermionic measurements operators factories
 */
#ifndef FERMION_MEAS_FACT_
#define FERMION_MEAS_FACT_

#include "include/pugi_interface.h"
#include "include/common_fields.hpp"
#include "Tools/RAIIFactory.hpp"

#include "quark_propagators.hpp"
#include "qprop.h"
#include "qprop_EvenOdd.h"
#include "qprop_DomainWall.hpp"
#include "Dirac_ops/dirac_Operator_Factory.hpp"
#include "Solver/solver_Factory.hpp"

/*!
 * @brief Abstract base class for creating QuarkPropagators
 */
class QuarkPropagatorFactory {
public:
  virtual QuarkPropagator* getQuarkProp(GaugeField& Field) = 0;
  virtual ~QuarkPropagatorFactory(){}
};

/*
class GenericMeasurementFactory {
public:
  virtual GenericMeasurement* getGMeas(GaugeField& Field) = 0;
  ~GenericMeasurementFactory();
}
*/

/////////////////////////////////////////////////
/*!
 @brief Concrete class for creating Quark Propagator Qprop operator
*/
class QPropFactory : public QuarkPropagatorFactory {
  RaiiFactoryObj<DiracWilsonLikeOperatorFactory> DiracObj;
  RaiiFactoryObj<SolverOperatorFactory> SolverObj;

  RaiiFactoryObj<DiracWilsonLike> Kernel;
  RaiiFactoryObj<Solver> Solv;
  const XML::node Qprop_node;

public:
  QPropFactory(XML::node node):Qprop_node(node){
    //some leaking expected - check!
    XML::descend(node,"Kernel");
    DiracObj.save(DiracOperators::createDiracWilsonLikeOperatorFactory(node));
    XML::next_sibling(node,"Solver");
    SolverObj.save(SolverOperators::createSolverOperatorFactory(node));
  }

  QuarkPropagator* getQuarkProp(GaugeField& Field){
    Kernel.save(DiracObj.get()->getDiracOperator(&(Field.U)));
    Solv.save(SolverObj.get()->getSolver(new Fopr_DdagD(Kernel.get())));

    return new Qprop(Kernel.get(), Solv.get());
  }
};

/*!
 @brief Concrete class for creating Quark Propagator Qprop_EvenOdd operator
*/
class QPropFactory_EvenOdd : public QuarkPropagatorFactory {
  RaiiFactoryObj<DiracWilsonEvenOddFactory> DiracObj;
  RaiiFactoryObj<SolverOperatorFactory> SolverObj;

  RaiiFactoryObj<DiracWilsonLike_EvenOdd> Kernel;
  RaiiFactoryObj<Solver> Solv;
  const XML::node Qprop_node;

public:
  QPropFactory_EvenOdd(XML::node node):Qprop_node(node){
    //some leaking expected - check!
    XML::descend(node,"Kernel");
    DiracObj.save(new DiracWilsonEvenOddFactory(node));
    XML::next_sibling(node,"Solver");
    SolverObj.save(SolverOperators::createSolverOperatorFactory(node));
  }

  QuarkPropagator* getQuarkProp(GaugeField& Field){
    Kernel.save(DiracObj.get()->getDiracOperator(&(Field.U)));

    std::cout<<"&DiracOp="<<Kernel.get() <<std::endl;

    Solv.save(SolverObj.get()->getSolver(new Fopr_DdagD(Kernel.get())));

    //std::cout<<"&Solv="<<Solv.get() <<std::endl;    
    //std::cout<<"Qprop being created "<<std::endl;
    return new Qprop_EvenOdd(Kernel.get(), Solv.get());
  }
};

/*!
 @brief Concrete class for creating Quark Propagator QpropDWF operator
*/
class QPropDWFFactory : public QuarkPropagatorFactory {
  RaiiFactoryObj<DiracDWF4dOperatorFactory> DiracFactory_;
  RaiiFactoryObj<Dirac_optimalDomainWall_4D> DWF4D_;
  const XML::node Qprop_node_;

public:
  QPropDWFFactory(XML::node node):Qprop_node_(node){
    XML::descend(node,"Kernel4d");
    DiracFactory_.save(DiracOperators::createDiracDWF4dOperatorFactory(node));
  }

  QuarkPropagator* getQuarkProp(GaugeField& Field){
    return getQuarkPropDW(Field);}

  QpropDWF* getQuarkPropDW(GaugeField& Field){
    DWF4D_.save(DiracFactory_.get()->getDiracOperator4D(&(Field.U)));
    return new QpropDWF(*DWF4D_.get());
  }
};

///////////////////////////////////////////////////

namespace QuarkPropagators{
  QuarkPropagatorFactory* createQuarkPropagatorFactory(const XML::node);
}


#endif 
