/*! @file dirac_Operator_Factory.hpp 
 *  @brief Declaration of Dirac operators factories
 Time-stamp: <2013-04-23 11:49:48 noaki>
 */
#ifndef DIRAC_FACT_
#define DIRAC_FACT_

#include "include/pugi_interface.h"
#include "Tools/RAIIFactory.hpp"
#include "dirac_wilson.hpp"
#include "dirac_wilson_EvenOdd.hpp"
#include "dirac_wilson_Brillouin.hpp"
#include "dirac_clover.hpp"
#include "dirac_DomainWall_4D_fullSolv.hpp"
#include "dirac_DomainWall_4D_eoSolv.hpp"
#include "dirac_DomainWall.hpp"
#include "dirac_DomainWall_EvenOdd.hpp"
#include "dirac_staggered_EvenOdd.hpp"
#include "dirac_staggered_EvenOdd_Adjoint.hpp"
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
    return getDiracOperatorWL(Gfield);}
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
  Dstagg::Dtype dt_;
public:
  DiracStaggeredEvenOddLikeOperatorFactory(Dstagg::Dtype dt):dt_(dt){}
  virtual DiracStaggeredEvenOddLike* 
  getDiracOperatorSTG(Dstagg::Dtype,Field* const) = 0;
  Dirac* getDiracOperator(Field* const Gfield){
    return getDiracOperatorSTG(dt_,Gfield);} // DdagDee is chosen as default
  virtual ~DiracStaggeredEvenOddLikeOperatorFactory(){}
};

//// conclete classes ////

/*! @brief Concrete class for creating Dirac Wilson operators */
class DiracWilsonFactory : public DiracWilsonLikeOperatorFactory {
  const XML::node Dirac_node_;
public:
  DiracWilsonFactory(const XML::node node):Dirac_node_(node){}
  Dirac_Wilson* getDiracOperatorWL(Field* const Gfield);
};

/*! @brief Concrete class for creating Dirac Wilson EvenOdd operators */
class DiracWilsonEvenOddFactory : public DiracWilsonLikeOperatorFactory {
  const XML::node Dirac_node_;
public:
  DiracWilsonEvenOddFactory(const XML::node node):Dirac_node_(node){}
  Dirac_Wilson_EvenOdd* getDiracOperatorWL(Field* const Gfield);
};

/*! @brief Concrete class for creating Dirac Wilson Brillouin operators */
class DiracWilsonBrillouinFactory : public DiracWilsonLikeOperatorFactory {
  const XML::node Dirac_node_;
  ImpType type_;
public:
  DiracWilsonBrillouinFactory(const XML::node node,ImpType type=Standard)
    :Dirac_node_(node),type_(type){}
  Dirac_Wilson_Brillouin* getDiracOperatorWL(Field* const Gfield);
};

/*! @brief Concrete class for creating Dirac Clover operators */
class DiracCloverFactory : public DiracWilsonLikeOperatorFactory {
  const XML::node Dirac_node_;
public:
  DiracCloverFactory(const XML::node node):Dirac_node_(node){}
  Dirac_Clover* getDiracOperatorWL(Field* const Gfield);
};

/*! @brief Concrete class for creating Dirac Optimal DWF-5d operators */
class DiracDWF5dFactory : public DiracDWF5dOperatorFactory {
  const XML::node Dirac_node_;
  RaiiFactoryObj<DiracWilsonLikeOperatorFactory> KernelFactory_;
public:
  DiracDWF5dFactory(XML::node node);

  Dirac_optimalDomainWall* getDiracOperatorWL(Field* const Gfield);
  Dirac_optimalDomainWall* getDiracOperatorPV(Field* const Gfield);
};

/*! @brief Concrete class for creating Dirac Optimal DWF-5d e/o operators */
class DiracDWF5dEvenOddFactory : public DiracDWF5dOperatorFactory {
  const XML::node Dirac_node_;
public:
  DiracDWF5dEvenOddFactory(XML::node node):Dirac_node_(node){}
  Dirac_optimalDomainWall_EvenOdd* getDiracOperatorWL(Field* const Gfield);
  Dirac_optimalDomainWall_EvenOdd* getDiracOperatorPV(Field* const Gfield);
};

/*! @brief Concrete class for creating Dirac_optimalDomainWall_4D_fullSolv */
class DiracDWF4DfullFactory : public DiracDWF4dOperatorFactory{
  const XML::node Dirac_node_;
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

public:
  DiracDWF4DfullFactory(XML::node node)
    :Dirac_node_(node){
    XML::descend(node,"Kernel5d", MANDATORY);
    DiracFactory_.save(new DiracDWF5dFactory(node));
    XML::next_sibling(node, "SolverDWF", MANDATORY);
    SolverFactory_.save(SolverOperators::createSolverOperatorFactory(node));
  }
  DiracWilsonLike* getDiracOperatorWL(Field* const Gfield);
  Dirac_optimalDomainWall_4D* getDiracOperator4D(Field* const Gfield);
};

/*! @brief Concrete class for creating Dirac_optimalDomainWall_4D_eoSolv*/
class DiracDWF4DeoFactory : public DiracDWF4dOperatorFactory{
  const XML::node Dirac_node_;
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

public:
  DiracDWF4DeoFactory(XML::node node);
  DiracWilsonLike* getDiracOperatorWL(Field* const Gfield);
  Dirac_optimalDomainWall_4D* getDiracOperator4D(Field* const Gfield);
};

/*! @brief Concrete class for creating Dirac_staggered_EvenOdd */
class DiracStaggeredEvenOddFactory
  : public DiracStaggeredEvenOddLikeOperatorFactory{
  const XML::node Dirac_node_;  
public:
  DiracStaggeredEvenOddFactory(const XML::node node,
			       Dstagg::Dtype dt=Dstagg::DdagDee)
    :Dirac_node_(node),DiracStaggeredEvenOddLikeOperatorFactory(dt){}
  Dirac_staggered_EvenOdd* getDiracOperatorSTG(Dstagg::Dtype dt,
					       Field* const Gfield);
};

/*! @brief Concrete class for creating Dirac_staggered_EvenOdd_Adjoint */
//Valid only in the case of NC=3
#if NC_==3
class DiracStaggeredEvenOddAdjointFactory
  : public DiracStaggeredEvenOddLikeOperatorFactory{
  const XML::node Dirac_node_;  
public:
  DiracStaggeredEvenOddAdjointFactory(const XML::node node,
				      Dstagg::Dtype dt=Dstagg::DdagDee)
    :Dirac_node_(node),DiracStaggeredEvenOddLikeOperatorFactory(dt){}
  Dirac_staggered_EvenOdd_Adjoint* getDiracOperatorSTG(Dstagg::Dtype dt,
						       Field* const Gfield);
};
#endif

//Add new factories here
//....

#endif
