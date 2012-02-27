/*!
 * @file dirac_Operator_Factory.hpp 
 *
 * @brief Declaration of Dirac operators factories
 */

#ifndef DIRAC_FACT_
#define DIRAC_FACT_

#include "include/pugi_interface.h"
#include "Tools/RAIIFactory.hpp"
#include "dirac_wilson.hpp"
#include "dirac_wilson_EvenOdd.hpp"
#include "dirac_clover.hpp"
#include "dirac_DomainWall_4D_fullSolv.hpp"
#include "dirac_DomainWall_4D_eoSolv.hpp"
#include "dirac_DomainWall.hpp"
#include "dirac_DomainWall_EvenOdd.hpp"
#include "Solver/solver_Factory.hpp"
#include "eoUtils.hpp"

/*!
 * @brief Abstract base class for creating Dirac operators
 */
class DiracWilsonLikeOperatorFactory {
public:
  virtual DiracWilsonLike* getDiracOperator(Field* const) = 0;
};

class DiracDWF5dOperatorFactory :public DiracWilsonLikeOperatorFactory{
public:
  virtual DiracWilsonLike* getDiracOperatorPV(Field* const) = 0;
};

class DiracDWF4dOperatorFactory :public DiracWilsonLikeOperatorFactory{
public:
  virtual Dirac_optimalDomainWall_4D* getDiracOperator4D(Field* const) = 0;
};
/////////////////////////////////////////////////////////
/*!
 * @brief Concrete class for creating Dirac Wilson operators
 */
class DiracWilsonFactory : public DiracWilsonLikeOperatorFactory {
  const XML::node Dirac_node;
public:
  DiracWilsonFactory(const XML::node node):Dirac_node(node){}

  Dirac_Wilson* getDiracOperator(Field* const GaugeField){
    return new Dirac_Wilson(Dirac_node,GaugeField);
  }
};
/////////////////////////////////////////////////////////
/*!
 * @brief Concrete class for creating Dirac Wilson EvenOdd operators
 */
class DiracWilsonEvenOddFactory : public DiracWilsonLikeOperatorFactory {
  const XML::node Dirac_node;
public:
  DiracWilsonEvenOddFactory(const XML::node node):Dirac_node(node){}

  Dirac_Wilson_EvenOdd* getDiracOperator(Field* const GaugeField){
    return new Dirac_Wilson_EvenOdd(Dirac_node,GaugeField);
  }
};
/////////////////////////////////////////////////////////
/*!
 * @brief Concrete class for creating Dirac Clover operators
 */
class DiracCloverFactory : public DiracWilsonLikeOperatorFactory {
  const XML::node Dirac_node;
public:
  DiracCloverFactory(const XML::node node):Dirac_node(node){}

  Dirac_Clover* getDiracOperator(Field* const GaugeField){
    return new Dirac_Clover(Dirac_node,GaugeField);
  }
};
/////////////////////////////////////////////////////////
/*!
 * @brief Concrete class for creating Dirac Optimal DWF-5d operators
 */
class DiracDWF5dFactory : public DiracDWF5dOperatorFactory {
  const XML::node Dirac_node;
public:
  DiracDWF5dFactory(XML::node node):Dirac_node(node){}

  Dirac_optimalDomainWall* getDiracOperator(Field* const GaugeField){
    return new Dirac_optimalDomainWall(Dirac_node,GaugeField);
  }
  Dirac_optimalDomainWall* getDiracOperatorPV(Field* const GaugeField){
    return new Dirac_optimalDomainWall(Dirac_node,GaugeField,PauliVillars);
  }
  ~DiracDWF5dFactory(){}
};
/////////////////////////////////////////////////////////
/*!
 * @brief Concrete class for creating Dirac Optimal DWF-5d e/o operators
 */
class DiracDWF5dEvenOddFactory : public DiracDWF5dOperatorFactory {
  const XML::node Dirac_node;
public:
  DiracDWF5dEvenOddFactory(XML::node node):Dirac_node(node){}

  Dirac_optimalDomainWall_EvenOdd* 
  getDiracOperator(Field* const GaugeField){
    return new Dirac_optimalDomainWall_EvenOdd(Dirac_node,GaugeField);
  }

  Dirac_optimalDomainWall_EvenOdd* 
  getDiracOperatorPV(Field* const GaugeField){
    return new 
      Dirac_optimalDomainWall_EvenOdd(Dirac_node,GaugeField,PauliVillars);
  }
  ~DiracDWF5dEvenOddFactory(){}
};
/////////////////////////////////////////////////////////
/*!
 * @brief Concrete class for creating Dirac_optimalDomainWall_4D_fullSolv
 */
class DiracDWF4DfullFactory : public DiracDWF4dOperatorFactory{
  // Factories
  RaiiFactoryObj<DiracDWF5dFactory> DiracFactory_;
  RaiiFactoryObj<SolverOperatorFactory> SolverFactory_;
  // Objects (Dodwf)
  RaiiFactoryObj<Dirac_optimalDomainWall> DW5D_;
  RaiiFactoryObj<Fopr_DdagD_Precondition> Fopr_;
  RaiiFactoryObj<Solver> Solver_;
  // Objects (PauliVillars)
  RaiiFactoryObj<Dirac_optimalDomainWall> DW5D_PV_;
  RaiiFactoryObj<Fopr_DdagD_Precondition> Fopr_PV_;
  RaiiFactoryObj<Solver> Solver_PV_;
  
  const XML::node Dirac_node_;

public:
  DiracDWF4DfullFactory(XML::node node)
    :Dirac_node_(node){
    XML::descend(node,"Kernel5d", MANDATORY);
    DiracFactory_.save(new DiracDWF5dFactory(node));
    XML::next_sibling(node, "SolverDWF", MANDATORY);
    SolverFactory_.save(SolverOperators::createSolverOperatorFactory(node));
  }

  DiracWilsonLike* getDiracOperator(Field* const GaugeField){
    return getDiracOperator4D(GaugeField);
  }

  Dirac_optimalDomainWall_4D* getDiracOperator4D(Field* const GaugeField){
    DW5D_.save(DiracFactory_.get()->getDiracOperator(GaugeField));
    Fopr_.save(new Fopr_DdagD_Precondition(DW5D_.get()));
    Solver_.save(SolverFactory_.get()->getSolver(Fopr_.get()));

    DW5D_PV_.save(DiracFactory_.get()->getDiracOperatorPV(GaugeField));
    Fopr_PV_.save(new Fopr_DdagD_Precondition(DW5D_PV_.get()));
    Solver_PV_.save(SolverFactory_.get()->getSolver(Fopr_PV_.get()));

    return new Dirac_optimalDomainWall_4D_fullSolv(DW5D_.get(),
						   DW5D_PV_.get(),
						   Solver_.get(),
						   Solver_PV_.get());
  }
  ~DiracDWF4DfullFactory(){}
};

/*!
 * @brief Concrete class for creating Dirac_optimalDomainWall_4D_eoSolv
 */
class DiracDWF4DeoFactory : public DiracDWF4dOperatorFactory{
  // Factories
  RaiiFactoryObj<DiracDWF5dFactory> DiracFactory_;
  RaiiFactoryObj<DiracDWF5dEvenOddFactory> DiracEOFactory_;
  RaiiFactoryObj<SolverOperatorFactory> SolverFactory_;
  // Objects (Dodwf)
  RaiiFactoryObj<Dirac_optimalDomainWall> DW5D_;
  RaiiFactoryObj<Dirac_optimalDomainWall_EvenOdd> DW5D_EO_;
  RaiiFactoryObj<Fopr_DdagD> FoprEO_;
  RaiiFactoryObj<Solver> SolverEO_;
  RaiiFactoryObj<EvenOddUtils::Inverter_WilsonLike> Inv_;
  // Objects (PauliVillars)
  RaiiFactoryObj<Dirac_optimalDomainWall> DW5DPV_;
  RaiiFactoryObj<Dirac_optimalDomainWall_EvenOdd> DW5D_EO_PV_;
  RaiiFactoryObj<Fopr_DdagD> FoprEO_PV_;
  RaiiFactoryObj<Solver> SolverEO_PV_;
  RaiiFactoryObj<EvenOddUtils::Inverter_WilsonLike> Inv_PV_;

  const XML::node Dirac_node_;

public:
  DiracDWF4DeoFactory(XML::node node)
    :Dirac_node_(node){
    XML::descend(node,"Kernel5d", MANDATORY);
    DiracFactory_.save(new DiracDWF5dFactory(node));
    DiracEOFactory_.save(new DiracDWF5dEvenOddFactory(node));
    XML::next_sibling(node, "SolverDWF", MANDATORY);
    SolverFactory_.save(SolverOperators::createSolverOperatorFactory(node));
  }

   DiracWilsonLike* getDiracOperator(Field* const GaugeField){
    return getDiracOperator4D(GaugeField);
  }

  Dirac_optimalDomainWall_4D* getDiracOperator4D(Field* const GaugeField){
    DW5D_.save(   DiracFactory_.get()->getDiracOperator(GaugeField));
    DW5D_EO_.save(DiracEOFactory_.get()->getDiracOperator(GaugeField));
    FoprEO_.save(new Fopr_DdagD(DW5D_EO_.get()));
    SolverEO_.save(SolverFactory_.get()->getSolver(FoprEO_.get()));
    Inv_.save(new EvenOddUtils::Inverter_WilsonLike(DW5D_EO_.get(),
						    SolverEO_.get()));

    DW5DPV_.save(    DiracFactory_.get()->getDiracOperatorPV(GaugeField));
    DW5D_EO_PV_.save(DiracEOFactory_.get()->getDiracOperatorPV(GaugeField));
    FoprEO_PV_.save(new Fopr_DdagD(DW5D_EO_PV_.get()));
    SolverEO_PV_.save(SolverFactory_.get()->getSolver(FoprEO_PV_.get()));
    Inv_PV_.save(new EvenOddUtils::Inverter_WilsonLike(DW5D_EO_PV_.get(),
						       SolverEO_PV_.get()));

    return new Dirac_optimalDomainWall_4D_eoSolv(DW5D_.get(),DW5DPV_.get(),
						 Inv_.get(),Inv_PV_.get());
  }
  ~DiracDWF4DeoFactory(){}
};

//Add new factories here
//....

//////////////////////////////////////////////////////////////
namespace DiracOperators {
  DiracWilsonLikeOperatorFactory* 
  createDiracWilsonLikeOperatorFactory(const XML::node);
  DiracDWF4dOperatorFactory* 
  createDiracDWF4dOperatorFactory(const XML::node);
  DiracDWF5dOperatorFactory* 
  createDiracDWF5dOperatorFactory(const XML::node);
}

#endif
