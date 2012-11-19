/*! @file dirac_Operator_Factory.hpp 
 *  @brief Declaration of Dirac operators factories
 */
#ifndef DIRAC_FACT_
#define DIRAC_FACT_

#include "include/pugi_interface.h"
#include "Tools/RAIIFactory.hpp"
#include "dirac_wilson.hpp"
#include "dirac_wilson_EvenOdd.hpp"
#include "dirac_wilson_Brillouin.hpp"
#include "dirac_wilson_Brillouin_Imp.hpp"
#include "dirac_clover.hpp"
#include "dirac_DomainWall_4D_fullSolv.hpp"
#include "dirac_DomainWall_4D_eoSolv.hpp"
#include "dirac_DomainWall.hpp"
#include "dirac_DomainWall_EvenOdd.hpp"
#include "dirac_staggered_EvenOdd.hpp"
#include "Solver/solver_Factory.hpp"
#include "eoUtils.hpp"

/*! @brief Abstract base class for creating Dirac operators */

class DiracOperatorFactory {
public:
  virtual Dirac* getDiracOperator(Field* const) = 0;
  virtual ~DiracOperatorFactory(){}
};

class DiracWilsonLikeOperatorFactory :public DiracOperatorFactory{
public:
  virtual DiracWilsonLike* getDiracOperatorWL(Field* const) = 0;
  Dirac* getDiracOperator(Field* const Gfield){
    return getDiracOperatorWL(Gfield);
  }
  virtual ~DiracWilsonLikeOperatorFactory(){}
};

class DiracDWF5dOperatorFactory :public DiracWilsonLikeOperatorFactory{
public:
  virtual DiracWilsonLike* getDiracOperatorPV(Field* const) = 0;
  virtual ~DiracDWF5dOperatorFactory(){}
};

class DiracDWF4dOperatorFactory :public DiracWilsonLikeOperatorFactory{
public:
  virtual Dirac_optimalDomainWall_4D* getDiracOperator4D(Field* const) = 0;
  virtual ~DiracDWF4dOperatorFactory(){}
};

class DiracStaggeredEvenOddLikeOperatorFactory :public DiracOperatorFactory{
public:
  virtual DiracStaggeredEvenOddLike* getDiracOperatorSTG(Field* const) = 0;
  Dirac* getDiracOperator(Field* const Gfield){
    return getDiracOperatorSTG(Gfield);
  }
  virtual ~DiracStaggeredEvenOddLikeOperatorFactory(){}
};


/*! @brief Concrete class for creating Dirac Wilson operators */
class DiracWilsonFactory : public DiracWilsonLikeOperatorFactory {
  const XML::node Dirac_node;
public:
  DiracWilsonFactory(const XML::node node):Dirac_node(node){}

  Dirac_Wilson* getDiracOperatorWL(Field* const Gfield){
    return new Dirac_Wilson(Dirac_node,Gfield);
  }
};

/*! @brief Concrete class for creating Dirac Wilson EvenOdd operators */
class DiracWilsonEvenOddFactory : public DiracWilsonLikeOperatorFactory {
  const XML::node Dirac_node;
public:
  DiracWilsonEvenOddFactory(const XML::node node):Dirac_node(node){}

  Dirac_Wilson_EvenOdd* getDiracOperatorWL(Field* const Gfield){
    return new Dirac_Wilson_EvenOdd(Dirac_node,Gfield);
  }
};

/*! @brief Concrete class for creating Dirac Wilson Brillouin operators */
class DiracWilsonBrillouinFactory : public DiracWilsonLikeOperatorFactory {
  const XML::node Dirac_node;
public:
  DiracWilsonBrillouinFactory(const XML::node node):Dirac_node(node){}

  Dirac_Wilson_Brillouin* getDiracOperatorWL(Field* const Gfield){
    return new Dirac_Wilson_Brillouin(Dirac_node,Gfield);
  }
};

/*! @brief Concrete class for creating Dirac Wilson Brillouin Improved operators */
class DiracWilsonBrillouinImpFactory : public DiracWilsonLikeOperatorFactory {
  const XML::node Dirac_node;
public:
  DiracWilsonBrillouinImpFactory(const XML::node node):Dirac_node(node){}

  Dirac_Wilson_Brillouin_Imp* getDiracOperatorWL(Field* const Gfield){
    return new Dirac_Wilson_Brillouin_Imp(Dirac_node,Gfield);
  }
};

/*! @brief Concrete class for creating Dirac Clover operators */
class DiracCloverFactory : public DiracWilsonLikeOperatorFactory {
  const XML::node Dirac_node;
public:
  DiracCloverFactory(const XML::node node):Dirac_node(node){}

  Dirac_Clover* getDiracOperatorWL(Field* const Gfield){
    return new Dirac_Clover(Dirac_node,Gfield);
  }
};

/*! @brief Concrete class for creating Dirac Optimal DWF-5d operators */
class DiracDWF5dFactory : public DiracDWF5dOperatorFactory {
  const XML::node Dirac_node;
public:
  DiracDWF5dFactory(XML::node node):Dirac_node(node){}

  Dirac_optimalDomainWall* getDiracOperatorWL(Field* const Gfield){
    return new Dirac_optimalDomainWall(Dirac_node,Gfield);
  }
  Dirac_optimalDomainWall* getDiracOperatorPV(Field* const Gfield){
    return new Dirac_optimalDomainWall(Dirac_node,Gfield,PauliVillars);
  }
  ~DiracDWF5dFactory(){}
};

/*! @brief Concrete class for creating Dirac Optimal DWF-5d e/o operators */
class DiracDWF5dEvenOddFactory : public DiracDWF5dOperatorFactory {
  const XML::node Dirac_node;
public:
  DiracDWF5dEvenOddFactory(XML::node node):Dirac_node(node){}

  Dirac_optimalDomainWall_EvenOdd* getDiracOperatorWL(Field* const Gfield){
    return new Dirac_optimalDomainWall_EvenOdd(Dirac_node,Gfield);
  }

  Dirac_optimalDomainWall_EvenOdd* 
  getDiracOperatorPV(Field* const Gfield){
    return new 
      Dirac_optimalDomainWall_EvenOdd(Dirac_node,Gfield,PauliVillars);
  }
  ~DiracDWF5dEvenOddFactory(){}
};

/*! @brief Concrete class for creating Dirac_optimalDomainWall_4D_fullSolv */
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

  DiracWilsonLike* getDiracOperatorWL(Field* const Gfield){
    return getDiracOperator4D(Gfield); }

  Dirac_optimalDomainWall_4D* getDiracOperator4D(Field* const Gfield){
    DW5D_.save(DiracFactory_.get()->getDiracOperatorWL(Gfield));
    Fopr_.save(new Fopr_DdagD_Precondition(DW5D_.get()));
    Solver_.save(SolverFactory_.get()->getSolver(Fopr_.get()));

    DW5D_PV_.save(DiracFactory_.get()->getDiracOperatorPV(Gfield));
    Fopr_PV_.save(new Fopr_DdagD_Precondition(DW5D_PV_.get()));
    Solver_PV_.save(SolverFactory_.get()->getSolver(Fopr_PV_.get()));

    return new Dirac_optimalDomainWall_4D_fullSolv(DW5D_.get(),
						   DW5D_PV_.get(),
						   Solver_.get(),
						   Solver_PV_.get());
  }
  ~DiracDWF4DfullFactory(){}
};

/*! @brief Concrete class for creating Dirac_optimalDomainWall_4D_eoSolv*/
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

   DiracWilsonLike* getDiracOperatorWL(Field* const Gfield){
    return getDiracOperator4D(Gfield); }

  Dirac_optimalDomainWall_4D* getDiracOperator4D(Field* const Gfield){
    DW5D_.save(   DiracFactory_.get()->getDiracOperatorWL(Gfield));
    DW5D_EO_.save(DiracEOFactory_.get()->getDiracOperatorWL(Gfield));
    FoprEO_.save(new Fopr_DdagD(DW5D_EO_.get()));
    SolverEO_.save(SolverFactory_.get()->getSolver(FoprEO_.get()));
    Inv_.save(new EvenOddUtils::Inverter_WilsonLike(DW5D_EO_.get(),
						    SolverEO_.get()));

    DW5DPV_.save(    DiracFactory_.get()->getDiracOperatorPV(Gfield));
    DW5D_EO_PV_.save(DiracEOFactory_.get()->getDiracOperatorPV(Gfield));
    FoprEO_PV_.save(new Fopr_DdagD(DW5D_EO_PV_.get()));
    SolverEO_PV_.save(SolverFactory_.get()->getSolver(FoprEO_PV_.get()));
    Inv_PV_.save(new EvenOddUtils::Inverter_WilsonLike(DW5D_EO_PV_.get(),
						       SolverEO_PV_.get()));

    return new Dirac_optimalDomainWall_4D_eoSolv(DW5D_.get(),DW5DPV_.get(),
						 Inv_.get(),Inv_PV_.get());
  }
  ~DiracDWF4DeoFactory(){}
};

/*! @brief Concrete class for creating Dirac_staggered_EvenOdd */
class DiracStaggeredEvenOddFactory
  : public DiracStaggeredEvenOddLikeOperatorFactory{

  const XML::node Dirac_node;
public:
  DiracStaggeredEvenOddFactory(const XML::node node):Dirac_node(node){}
  Dirac_staggered_EvenOdd* getDiracOperatorSTG(Field* const Gfield){
    return new Dirac_staggered_EvenOdd(Dirac_node,Gfield);
  }
};

//Add new factories here
//....

//////////////////////////////////////////////////////////////
namespace DiracOperators {
  DiracOperatorFactory* 
  createGeneralDiracOperatorFactory(const XML::node);

  DiracWilsonLikeOperatorFactory* 
  createDiracWilsonLikeOperatorFactory(const XML::node);

  DiracDWF4dOperatorFactory* 
  createDiracDWF4dOperatorFactory(const XML::node);

  DiracDWF5dOperatorFactory* 
  createDiracDWF5dOperatorFactory(const XML::node);

  DiracWilsonLikeOperatorFactory* 
  createGeneralDiracWilsonLikeOperatorFactory(const XML::node);
  
  DiracStaggeredEvenOddLikeOperatorFactory* 
  createDiracStaggeredEvenOddLikeOperatorFactory(const XML::node);
}

#endif
