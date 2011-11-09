/*!
 * @file fermion_meas_factory.hpp 
 *
 * @brief Fermionic measurements operators factories
 */

#ifndef FERMION_MEAS_FACT_
#define FERMION_MEAS_FACT_

#include "include/pugi_interface.h"
#include "include/common_fields.hpp"
#include "Tools/RAIIFactory.hpp"

#include "quark_propagators.hpp"
#include "qprop.h"
#include "qprop_optimalDomainWall.hpp"
#include "Dirac_ops/dirac_Factory.hpp"
#include "Solver/solver_Factory.hpp"


/*!
 * @brief Abstract base class for creating QuarkPropagators
 *
 */
class QuarkPropagatorCreator {
public:
  virtual QuarkPropagator* getQuarkProp(GaugeField& Field) = 0;
};

/////////////////////////////////////////////////
/*!
 @brief Concrete class for creating Quark Propagator Qprop operator
*/
class QPropCreator : public QuarkPropagatorCreator {
  RaiiFactoryObj<DiracOperatorCreator> DiracObj;
  RaiiFactoryObj<SolverOperatorCreator> SolverObj;

  RaiiFactoryObj<Dirac> Kernel;
  RaiiFactoryObj<Solver> Solv;
  const XML::node Qprop_node;

public:
  QPropCreator(XML::node node):
    Qprop_node(node)
    {
      //some leaking expected - check!
      XML::descend(node,"Kernel");
      DiracObj.save(DiracOperators::createDiracOperatorFactory(node));
      XML::next_sibling(node,"Solver");
      SolverObj.save(SolverOperators::createSolverOperatorFactory(node));
    };

  ~QPropCreator(){
  }

  QuarkPropagator* getQuarkProp(GaugeField& Field){
    Kernel.save(DiracObj.get()->getDiracOperator(&(Field.U)));
    Solv.save(SolverObj.get()->getSolver(new Fopr_DdagD(Kernel.get())));

    return new Qprop(Kernel.get(), Solv.get());
  };

};

/*!
 @brief Concrete class for creating Quark Propagator QpropDWF operator
*/
class QPropDWFCreator : public QuarkPropagatorCreator {
  RaiiFactoryObj<DiracDWF4dCreator> DiracObj; //very specific
  
  RaiiFactoryObj<Dirac_optimalDomainWall_4D> Kernel_DWF4D;
  const XML::node Qprop_node;

public:
  QPropDWFCreator(XML::node node):
    Qprop_node(node)
    {
      XML::descend(node,"Kernel4d");
      DiracObj.save(new DiracDWF4dCreator(node));
    };

  ~QPropDWFCreator(){}

  QuarkPropagator* getQuarkProp(GaugeField& Field){
    Kernel_DWF4D.save(DiracObj.get()->getDiracOperator(&(Field.U)));

    return new QpropDWF(*Kernel_DWF4D.get());
  };

};

///////////////////////////////////////////////////

namespace QuarkPropagators{
  QuarkPropagatorCreator* createQuarkPropagatorFactory(const XML::node);
}


#endif 
