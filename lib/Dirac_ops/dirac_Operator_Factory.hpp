/*! @file dirac_Operator_Factory.hpp 
 *  @brief Declaration of Dirac operators factories
 Time-stamp: <2013-05-23 10:28:58 noaki>
 */
#ifndef DIRAC_FACT_
#define DIRAC_FACT_

#include "include/pugi_interface.h"
#include "Tools/RAIIFactory.hpp"

#include "dirac_wilson_EvenOdd.hpp"
#include "dirac_wilson_Brillouin.hpp"
#include "dirac_clover.hpp"
#include "dirac_Mobius.hpp"
#include "dirac_DomainWall_4D_fullSolv.hpp"
#include "dirac_DomainWall_4D_eoSolv.hpp"
#include "dirac_DomainWall.hpp"
#include "dirac_DomainWall_EvenOdd.hpp"
#include "dirac_staggered_EvenOdd.hpp"
#include "dirac_staggered_EvenOdd_Adjoint.hpp"
class SolverFactory;

/*! @brief Abstract base class for creating Dirac operators */

class DiracFactory {
public:
  virtual Dirac* getDirac(Field* const) = 0;
  virtual ~DiracFactory(){}
};

class DiracWilsonLikeFactory :public DiracFactory{
public:
  virtual DiracWilsonLike* getDiracWL(Field* const) = 0;
  Dirac* getDirac(Field* const Gfield){
    return getDiracWL(Gfield);}
  virtual ~DiracWilsonLikeFactory(){}
};

class DiracWilsonLikeEvenOddFactory
  : public DiracWilsonLikeFactory{
public:
  /*!@brief this virtual function is to be wrapped and to give eo-/oe- operators to the higher level operator */
  virtual DiracWilsonLike_EvenOdd* getDiracEO(Field* const) = 0;

  DiracWilsonLike* getDiracWL(Field* const Gfield){
    return getDiracEO(Gfield);}
  virtual ~DiracWilsonLikeEvenOddFactory(){}
};

class DiracDWF5dFactory :public DiracWilsonLikeFactory{
public:
  virtual DiracWilsonLike* getDiracPV(Field* const) = 0;
  virtual ~DiracDWF5dFactory(){}
};

class DiracDWF4dFactory :public DiracWilsonLikeFactory{
public:
  virtual Dirac_optimalDomainWall_4D* getDirac4D(Field* const) = 0;
  DiracWilsonLike* getDiracWL(Field* const Gfield){
    return getDirac4D(Gfield);}

  virtual ~DiracDWF4dFactory(){}
};

class DiracStaggeredEvenOddLikeFactory :public DiracFactory{
  Dstagg::Dtype dt_;
public:
  DiracStaggeredEvenOddLikeFactory(Dstagg::Dtype dt):dt_(dt){}
  virtual DiracStaggeredEvenOddLike* 
  getDiracSTG(Dstagg::Dtype,Field* const) = 0;
  Dirac* getDirac(Field* const Gfield){
    return getDiracSTG(dt_,Gfield);} // DdagDee is chosen as default
  virtual ~DiracStaggeredEvenOddLikeFactory(){}
};

/////////////////////// conclete classes ///////////////////////////

/*! @brief Concrete class for creating Dirac Wilson operators */
class DiracWilsonFactory : public DiracWilsonLikeFactory {
  const XML::node Dirac_node_;
public:
  DiracWilsonFactory(const XML::node node):Dirac_node_(node){}
  Dirac_Wilson* getDiracWL(Field* const Gfield);
};

/*! @brief Concrete class for creating Dirac Wilson EvenOdd operators */
class DiracWilsonEvenOddFactory : public DiracWilsonLikeEvenOddFactory {
  const XML::node Dirac_node_;
public:
  DiracWilsonEvenOddFactory(const XML::node node):Dirac_node_(node){}
  Dirac_Wilson_EvenOdd* getDiracEO(Field* const Gfield);
};

/*! @brief Concrete class for creating Dirac Wilson Brillouin operators */
class DiracWilsonBrillouinFactory : public DiracWilsonLikeFactory {
  const XML::node Dirac_node_;
  ImpType type_;
public:
  DiracWilsonBrillouinFactory(const XML::node node,ImpType type=Standard)
    :Dirac_node_(node),type_(type){}
  Dirac_Wilson_Brillouin* getDiracWL(Field* const Gfield);
};

/*! @brief Concrete class for creating Dirac Clover operators */
class DiracCloverFactory : public DiracWilsonLikeFactory {
  const XML::node Dirac_node_;
public:
  DiracCloverFactory(const XML::node node):Dirac_node_(node){}
  Dirac_Clover* getDiracWL(Field* const Gfield);
};

/*! @brief Concrete class for creating Dirac Mobius operators */
class DiracMobiusFactory : public DiracWilsonLikeFactory{
  // Factories
  RaiiFactoryObj<DiracWilsonLikeFactory> DiracFactory_;
  RaiiFactoryObj<SolverFactory> SolverFactory_;
  // Objects (Dodwf)
  RaiiFactoryObj<DiracWilsonLike> D_;
  RaiiFactoryObj<Fopr_DdagD_Precondition> Fopr_;
  RaiiFactoryObj<Solver> Solver_;

  const XML::node Dirac_node_;
public:
  DiracMobiusFactory(XML::node node);
  Dirac_Mobius* getDiracWL(Field* const Gfield);
};

/*! @brief Concrete class for creating Dirac Optimal DWF-5d operators */
class DiracDomainWall5dFactory : public DiracDWF5dFactory {
  const XML::node Dirac_node_;
  RaiiFactoryObj<DiracWilsonLikeFactory> KernelFactory_;
public:
  DiracDomainWall5dFactory(XML::node node);

  Dirac_optimalDomainWall* getDiracWL(Field* const Gfield);
  Dirac_optimalDomainWall* getDiracPV(Field* const Gfield);
};

/*! @brief Concrete class for creating Dirac Optimal DWF-5d e/o operators */
class DiracDomainWall5dEvenOddFactory : public DiracDWF5dFactory {
  const XML::node Dirac_node_;
  RaiiFactoryObj<DiracWilsonLikeEvenOddFactory> KernelFactory_;
  RaiiFactoryObj<DiracWilsonLike_EvenOdd> Kernel_;
public:
  DiracDomainWall5dEvenOddFactory(XML::node node);
  Dirac_optimalDomainWall_EvenOdd* getDiracWL(Field* const Gfield);
  Dirac_optimalDomainWall_EvenOdd* getDiracPV(Field* const Gfield);
};

/*! @brief Concrete class for creating Dirac_optimalDomainWall_4D_fullSolv */
class DiracDWF4DfullFactory : public DiracDWF4dFactory{
  const XML::node Dirac_node_;
  // Factories
  RaiiFactoryObj<DiracDomainWall5dFactory> DiracFactory_;
  RaiiFactoryObj<SolverFactory> SolverFactory_;
  // Objects (Dodwf)
  RaiiFactoryObj<Dirac_optimalDomainWall> DW5D_;
  RaiiFactoryObj<Fopr_DdagD_Precondition> Fopr_;
  RaiiFactoryObj<Solver> Solver_;
  // Objects (PauliVillars)
  RaiiFactoryObj<Dirac_optimalDomainWall> DW5D_PV_;
  RaiiFactoryObj<Fopr_DdagD_Precondition> Fopr_PV_;
  RaiiFactoryObj<Solver> Solver_PV_;

public:
  DiracDWF4DfullFactory(XML::node node);
  Dirac_optimalDomainWall_4D* getDirac4D(Field* const Gfield);
};

/*! @brief Concrete class for creating Dirac_optimalDomainWall_4D_eoSolv*/
class DiracDWF4DeoFactory : public DiracDWF4dFactory{
  const XML::node Dirac_node_;
  // Factories
  RaiiFactoryObj<DiracDomainWall5dEvenOddFactory> DiracEOFactory_;
  RaiiFactoryObj<SolverFactory> SolverFactory_;
  // Objects (Dodwf)
  RaiiFactoryObj<Dirac_optimalDomainWall_EvenOdd> DW5D_EO_;
  RaiiFactoryObj<Fopr_DdagD> FoprEO_;
  RaiiFactoryObj<Solver> SolverEO_;
  RaiiFactoryObj<EvenOddUtils::Inverter_WilsonLike> Inv_;
  // Objects (PauliVillars)
  RaiiFactoryObj<Dirac_optimalDomainWall_EvenOdd> DW5D_EO_PV_;
  RaiiFactoryObj<Fopr_DdagD> FoprEO_PV_;
  RaiiFactoryObj<Solver> SolverEO_PV_;
  RaiiFactoryObj<EvenOddUtils::Inverter_WilsonLike> Inv_PV_;

public:
  DiracDWF4DeoFactory(XML::node node);
  Dirac_optimalDomainWall_4D* getDirac4D(Field* const Gfield);
};

/*! @brief Concrete class for creating Dirac_staggered_EvenOdd */
class DiracStaggeredEvenOddFactory
  : public DiracStaggeredEvenOddLikeFactory{
  const XML::node Dirac_node_;  
public:
  DiracStaggeredEvenOddFactory(const XML::node node,
			       Dstagg::Dtype dt=Dstagg::DdagDee)
    :Dirac_node_(node),DiracStaggeredEvenOddLikeFactory(dt){}
  Dirac_staggered_EvenOdd* getDiracSTG(Dstagg::Dtype dt,
					       Field* const Gfield);
};

/*! @brief Concrete class for creating Dirac_staggered_EvenOdd_Adjoint */
//Valid only in the case of NC=3
#if NC_==3
class DiracStaggeredEvenOddAdjointFactory
  : public DiracStaggeredEvenOddLikeFactory{
  const XML::node Dirac_node_;  
public:
  DiracStaggeredEvenOddAdjointFactory(const XML::node node,
				      Dstagg::Dtype dt=Dstagg::DdagDee)
    :Dirac_node_(node),DiracStaggeredEvenOddLikeFactory(dt){}
  Dirac_staggered_EvenOdd_Adjoint* getDiracSTG(Dstagg::Dtype dt,
						       Field* const Gfield);
};
#endif

//Add new factories here
//....

#endif
