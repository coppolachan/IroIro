/*!@file quark_prop_meas_factory.hpp 
  @brief QuarkPropagators measurements operators factories
*/
#ifndef QUARKPROP_MEAS_FACT_
#define QUARKPROP_MEAS_FACT_

#include "include/pugi_interface.h"
#include "Tools/RAIIFactory.hpp"
#include "inputConfig.hpp"
// Concrete objects definitions
#include "qprop.hpp"
#include "qprop_EvenOdd.hpp"
#include "qprop_DomainWall.hpp"
#include "Dirac_ops/dirac_Operator_Factory.hpp"

/*! @brief Abstract base class for creating QuarkPropagators */
class QuarkPropagatorFactory {
  virtual QuarkPropagator* createQuarkProp(InputConfig&)= 0;
public:
  QuarkPropagator* getQuarkProp(InputConfig& input){
    return createQuarkProp(input); }
  virtual ~QuarkPropagatorFactory(){}
};

/////////////////////////////////////////////////
/*! @brief Concrete class for creating Quark Propagator Qprop operator */

class QPropFactory : public QuarkPropagatorFactory {
  RaiiFactoryObj<DiracWilsonLikeFactory> Dfactory_;
  RaiiFactoryObj<SolverFactory> slvFactory_;

  RaiiFactoryObj<DiracWilsonLike> Kernel_;
  RaiiFactoryObj<Solver> Solv_;
  RaiiFactoryObj<Fopr_DdagD> DdagD_;
  
  QuarkPropagator* createQuarkProp(InputConfig&);
  XML::node node_;
public:
  QPropFactory(const XML::node&);
  QPropFactory(const XML::node&,double mass);
};

/*! @brief Concrete class for creating Quark Propagator Qprop_EvenOdd operator*/
class QPropFactory_EvenOdd : public QuarkPropagatorFactory {
  RaiiFactoryObj<DiracWilsonLikeEvenOddFactory> Dfactory_;
  RaiiFactoryObj<SolverFactory> slvFactory_;

  RaiiFactoryObj<DiracWilsonLike_EvenOdd> Kernel_;
  RaiiFactoryObj<Solver> Solv_;
  RaiiFactoryObj<Fopr_DdagD> DdagD_;

  QuarkPropagator* createQuarkProp(InputConfig&);
  XML::node node_;
public:
  QPropFactory_EvenOdd(const XML::node&);
  QPropFactory_EvenOdd(const XML::node&,double mass);
};

/*! @brief Concrete class for creating Quark Propagator QpropDWF operator */
class QPropDWFFactory : public QuarkPropagatorFactory {
  RaiiFactoryObj<DiracDWF4dFactory> Dfactory_;
  RaiiFactoryObj<Dirac_DomainWall_4D> DWF4D_;

  QpropDWF* createQuarkProp(InputConfig&);
  XML::node node_;
public:
  QPropDWFFactory(const XML::node&);
  QPropDWFFactory(const XML::node&,double mass);

  QpropDWF* getQuarkPropDW(InputConfig& input){createQuarkProp(input);}
};

namespace QuarkPropagators{
  QuarkPropagatorFactory* createQuarkPropagatorFactory(const XML::node&);
  QuarkPropagatorFactory* createQuarkPropagatorFactory(const XML::node&,
						       double mass);
}

#endif 
