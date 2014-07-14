/*! @file dirac_Operator_Factory.hpp 
 *  @brief Declaration of Dirac operators factories
 Time-stamp: <2014-07-04 17:49:11 noaki>
 */
#ifndef DIRAC_FACT_
#define DIRAC_FACT_

#include "pugi_interface.h"
#include "Tools/RAIIFactory.hpp"

#include "inputConfig.hpp"
#include "dirac_wilson_adjoint.hpp"
#include "dirac_wilson_EvenOdd.hpp"
#include "dirac_wilson_Brillouin.hpp"
#include "dirac_wilson_Brillouin_OSS.hpp"
#include "dirac_wilson_FiniteDensity.hpp"
#include "dirac_clover.hpp"
#include "dirac_Mobius.hpp"
#include "dirac_DomainWall_adjoint.hpp"
#include "dirac_DomainWall_4D_fullSolv.hpp"
#include "dirac_DomainWall_4D_eoSolv.hpp"
#include "dirac_DomainWall_EvenOdd.hpp"
#include "dirac_DomainWall_adjoint_EvenOdd.hpp"
#include "dirac_DWoverlap.hpp"
#include "dirac_LowModeDeflation_ExactEigen.hpp"
#include "dirac_LowModeDeflation_Approx.hpp"
#include "dirac_staggered_EvenOdd.hpp"
#include "dirac_staggered_EvenOdd_Adjoint.hpp"
#include "EigenModes/eigenModes.hpp"

class SolverFactory;
class SolverCG_DWF_opt_Factory;
/* employed NVI (non-virtual-interface) idiom: Effective C++ (item 35) */

/*! @brief Abstract base class for creating Dirac operators */
class DiracFactory {
public:
  virtual Dirac* getDirac(InputConfig&) = 0;
  virtual ~DiracFactory(){}
};

class DiracWilsonLikeFactory :public DiracFactory{
  virtual DiracWilsonLike* createDirac(InputConfig&) = 0;
public:
  DiracWilsonLike* getDirac(InputConfig& input){return createDirac(input);}
  virtual ~DiracWilsonLikeFactory(){}
};

class DiracWilsonLikeEvenOddFactory: virtual public DiracWilsonLikeFactory{
  virtual DiracWilsonLike_EvenOdd* createDirac(InputConfig&) = 0;
public:
  DiracWilsonLike_EvenOdd* getDirac(InputConfig& input){return createDirac(input);}
  virtual ~DiracWilsonLikeEvenOddFactory(){}
};

/*!@brief this abstruct class is used for HMC, not for measurements */
class DiracDWF5dFactory :virtual public DiracWilsonLikeFactory{
  virtual DiracWilsonLike* createDirac(InputConfig&) = 0;
  virtual DiracWilsonLike* createDiracPV(InputConfig&) =0;
public:
  DiracWilsonLike* getDirac(InputConfig& input){return createDirac(input);} 
  DiracWilsonLike* getDiracPV(InputConfig& input){return createDiracPV(input);}
  virtual ~DiracDWF5dFactory(){}
};

/*!@brief this abstruct class is used only for measurements */
class DiracDWF4dFactory :public DiracWilsonLikeFactory{
  virtual Dirac_DomainWall_4D* createDirac(InputConfig&) = 0;
public:
  Dirac_DomainWall_4D* getDirac(InputConfig& input){
    return createDirac(input);} //name surpression
  virtual ~DiracDWF4dFactory(){}
};

/*!@brief this abstruct class is used only for measurements */
class DiracDeflationFactory :public DiracWilsonLikeFactory{
  virtual Dirac_LowModeDeflation* createDirac(InputConfig&) = 0;
public:
  Dirac_LowModeDeflation* getDirac(InputConfig& input){return createDirac(input);}
  virtual ~DiracDeflationFactory(){}
};

class DiracStaggeredEvenOddLikeFactory :public DiracFactory{
  virtual DiracStaggeredEvenOddLike* createDirac(InputConfig&) = 0;
  virtual DiracStaggeredEvenOddLike* createDoo(InputConfig&) = 0;
public:
  DiracStaggeredEvenOddLike* getDirac(InputConfig& input){
    return createDirac(input);}//name surpression, DdagDee is chosen as default

  DiracStaggeredEvenOddLike* getDoo(InputConfig& input){
    return createDoo(input);}//name surpression, DdagDee is chosen as default
  
  virtual ~DiracStaggeredEvenOddLikeFactory(){}
};

/////////////////////// conclete classes ///////////////////////////
/*! @brief Concrete class for creating Dirac Wilson operators */
class DiracWilsonFactory : public DiracWilsonLikeFactory {
  const XML::node Dirac_node_;
  DiracWilsonLike* createDirac(InputConfig&);
public:
  DiracWilsonFactory(const XML::node& node):Dirac_node_(node){}
};

/////////////
/*! @brief Concrete class for creating Dirac Wilson (adjoint) operators */
class DiracWilsonAdjointFactory : public DiracWilsonLikeFactory {
  const XML::node Dirac_node_;
  DiracWilsonLike* createDirac(InputConfig&);
public:
  DiracWilsonAdjointFactory(const XML::node& node):Dirac_node_(node){}
};

/////////////
/*! @brief Concrete class for creating Dirac Wilson Brillouin operators */
class DiracWilsonBrillouinFactory : public DiracWilsonLikeFactory {
  XML::node Dirac_node_;
  ImpType type_;
  DiracWilsonLike* createDirac(InputConfig&);
public:
  DiracWilsonBrillouinFactory(const XML::node& node,ImpType type=Standard)
    :Dirac_node_(node),type_(type){}
};

/////////////
/*! @brief Concrete class for creating Dirac Wilson Brillouin OSS operators */
class DiracWilsonBrillouinOSSFactory : public DiracWilsonLikeFactory {
  XML::node Dirac_node_;
  ImpType type_;
  DiracWilsonLike* createDirac(InputConfig&);
public:
  DiracWilsonBrillouinOSSFactory(const XML::node& node,ImpType type=Standard)
    :Dirac_node_(node),type_(type){}
};
////////////
/*! @brief Concrete class for creating Dirac Wilson FiniteDensity operators */
class DiracWilsonFiniteDensityFactory : public DiracWilsonLikeFactory {
  const XML::node Dirac_node_;
  DiracWilsonLike* createDirac(InputConfig&);
public:
  DiracWilsonFiniteDensityFactory(const XML::node& node):Dirac_node_(node){}
};

/////////////
/*! @brief Concrete class for creating Dirac Clover operators */
class DiracCloverFactory : public DiracWilsonLikeFactory {

  RaiiFactoryObj<DiracWilsonLikeFactory> DiracFactory_;
  RaiiFactoryObj<DiracWilsonLike> D_;
  XML::node Dirac_node_;
  DiracWilsonLike* createDirac(InputConfig&);
public:
  DiracCloverFactory(XML::node node);
};

////////////
/*! @brief Concrete class for creating Dirac Mobius operators */
class DiracMobiusFactory : public DiracWilsonLikeFactory{
  // Factories
  RaiiFactoryObj<DiracWilsonLikeFactory> DiracFactory_;
  RaiiFactoryObj<SolverFactory> SolverFactory_;
  // Objects (Dodwf)
  RaiiFactoryObj<DiracWilsonLike> D_;
  RaiiFactoryObj<Fopr_DdagD> Fopr_;
  RaiiFactoryObj<Solver> Solver_;

  XML::node Dirac_node_;
  DiracWilsonLike* createDirac(InputConfig&);
public:
  DiracMobiusFactory(XML::node node);
};

/////////////
/*! @brief Concrete class for creating Dirac DWoverap operators */
class DiracDWoverlapFactory : public DiracWilsonLikeFactory {
  XML::node Dirac_node_;
  RaiiFactoryObj<DiracDWF4dFactory> DW4dFactory_;
  RaiiFactoryObj<Dirac_DomainWall_4D> DW4d_;
  
  DiracWilsonLike* createDirac(InputConfig&);
public:
  DiracDWoverlapFactory(XML::node node);
};

/////////////
/*! @brief Concrete class for creating Dirac_Deflation_ExactEigen */
class DiracDeflationExactFactory : public DiracDeflationFactory {
  XML::node Dirac_node_;
  RaiiFactoryObj<DiracWilsonLikeFactory> DwFactory_;
  RaiiFactoryObj<DiracWilsonLike> Dw_;

  Dirac_LowModeDeflation* createDirac(InputConfig&);
public:
  DiracDeflationExactFactory(XML::node node);
};

/////////////
/*! @brief Concrete class for creating Dirac_Deflation_Approx */
class DiracDeflationApproxFactory : public DiracDeflationFactory {
  XML::node Dirac_node_;
  RaiiFactoryObj<DiracWilsonLikeFactory> DwFactory_;
  RaiiFactoryObj<DiracWilsonLike> Dw_;

  Dirac_LowModeDeflation* createDirac(InputConfig&);
public:
  DiracDeflationApproxFactory(XML::node node);
};

/////////////
/*! @brief Concrete class for creating Dirac_Deflation_ApproxSubspace */
class DiracDeflationApproxSubspaceFactory : public DiracWilsonLikeFactory {
  XML::node Dirac_node_;
  Dirac_LowModeDeflation* createDirac(InputConfig&);
public:
  DiracDeflationApproxSubspaceFactory(const XML::node& node):Dirac_node_(node){}
};

/////////////
/*! @brief Concrete class for creating Dirac DWF-5d operators */
class DiracDomainWall5dFactory : public DiracDWF5dFactory {
  XML::node Dirac_node_;
  RaiiFactoryObj<DiracWilsonLikeFactory> KernelFactory_;
  RaiiFactoryObj<DiracWilsonLike> Kernel_;

  Dirac_DomainWall* createDirac(InputConfig&);
  Dirac_DomainWall* createDiracPV(InputConfig&);
public:
  DiracDomainWall5dFactory(XML::node node);
  Dirac_DomainWall* getDirac(InputConfig& input){
    return createDirac(input);}  //name surpression
  Dirac_DomainWall* getDiracPV(InputConfig& input){
    return createDiracPV(input);}  //name surpression
};

/*! @brief Concrete class for creating Dirac DWF_Adjoint-5d operators */
class DiracDomainWall5dAdjointFactory : public DiracDWF5dFactory {
  XML::node Dirac_node_;
  RaiiFactoryObj<DiracWilsonLikeFactory> KernelFactory_;
  RaiiFactoryObj<DiracWilsonLike> Kernel_;

  Dirac_DomainWall_Adjoint* createDirac(InputConfig&);
  Dirac_DomainWall_Adjoint* createDiracPV(InputConfig&);
public:
  DiracDomainWall5dAdjointFactory(XML::node node);
  Dirac_DomainWall_Adjoint* getDirac(InputConfig& input){
    return createDirac(input);}  //name surpression
  Dirac_DomainWall_Adjoint* getDiracPV(InputConfig& input){
    return createDiracPV(input);}  //name surpression
};

/////////////
/*! @brief Concrete class for creating Dirac Wilson EvenOdd operators */
class DiracWilsonEvenOddFactory : public DiracWilsonLikeEvenOddFactory {
  XML::node Dirac_node_;
  DiracWilsonLike_EvenOdd* createDirac(InputConfig&);
public:
  DiracWilsonEvenOddFactory(const XML::node& node):Dirac_node_(node){}
};

/////////////
/*! @brief Concrete class for creating Dirac Wilson Adjoint EvenOdd operators */
class DiracWilsonAdjointEvenOddFactory : public DiracWilsonLikeEvenOddFactory {
  XML::node Dirac_node_;
  DiracWilsonLike_EvenOdd* createDirac(InputConfig&);
public:
  DiracWilsonAdjointEvenOddFactory(const XML::node& node):Dirac_node_(node){}
};

//////////////
/*! @brief Concrete class for creating DWF-5d e/o operator as EvenOdd */
class DiracEvenOdd_DWF5dFactory : public DiracWilsonLikeEvenOddFactory,
				  public DiracDWF5dFactory{
  XML::node Dirac_node_;
  RaiiFactoryObj<DiracWilsonLikeEvenOddFactory> KernelFactory_;
  RaiiFactoryObj<DiracWilsonLike_EvenOdd> Kernel_;
  
  Dirac_DomainWall_EvenOdd* createDirac(InputConfig&);
  Dirac_DomainWall_EvenOdd* createDiracPV(InputConfig&);
public:
  DiracEvenOdd_DWF5dFactory(XML::node node);
  Dirac_DomainWall_EvenOdd* getDirac(InputConfig& input){
    return createDirac(input);}  //name surpression
  Dirac_DomainWall_EvenOdd* getDiracPV(InputConfig& input){
    return createDiracPV(input);}  //name surpression
};

//////////////
/*! @brief Concrete class for creating DWFadjoint-5d e/o operator as EvenOdd */
class DiracEvenOdd_DWF5dAdjointFactory : public DiracWilsonLikeEvenOddFactory,
					 public DiracDWF5dFactory{
  XML::node Dirac_node_;
  RaiiFactoryObj<DiracWilsonLikeEvenOddFactory> KernelFactory_;
  RaiiFactoryObj<DiracWilsonLike_EvenOdd> Kernel_;
  
  Dirac_DomainWall_Adjoint_EvenOdd* createDirac(InputConfig&);
  Dirac_DomainWall_Adjoint_EvenOdd* createDiracPV(InputConfig&);
public:
  DiracEvenOdd_DWF5dAdjointFactory(XML::node node);
  Dirac_DomainWall_Adjoint_EvenOdd* getDirac(InputConfig& input){
    return createDirac(input);}   //name surpression
  Dirac_DomainWall_Adjoint_EvenOdd* getDiracPV(InputConfig& input){
    return createDiracPV(input);} //name surpression
};

//////////////
/*! @brief Concrete class for creating Dirac_DomainWall_4D_fullSolv */
class DiracDWF4DfullFactory : public DiracDWF4dFactory{
  XML::node Dirac_node_;
  // Factories
  RaiiFactoryObj<DiracDomainWall5dFactory> DiracFactory_;
  RaiiFactoryObj<SolverFactory> SolverFactory_;
  // Objects (Dodwf)
  RaiiFactoryObj<Dirac_DomainWall> DW5D_;
  RaiiFactoryObj<Fopr_DdagD> Fopr_;
  RaiiFactoryObj<Solver> Solver_;
  // Objects (PauliVillars)
  RaiiFactoryObj<Dirac_DomainWall> DW5dPV_;
  RaiiFactoryObj<Fopr_DdagD> FoprPV_;
  RaiiFactoryObj<Solver> SolverPV_;

  Dirac_DomainWall_4D* createDirac(InputConfig&);
  
  DW5dPrecond prec_;
public:
  DiracDWF4DfullFactory(XML::node node,DW5dPrecond prec=NoPrecond);
};

//////////////
/*! @brief Concrete class for creating Dirac_DomainWall_4D_eoSolv*/
class DiracDWF4DeoFactory : public DiracDWF4dFactory{
  XML::node Dirac_node_;
  // Factories
  RaiiFactoryObj<DiracEvenOdd_DWF5dFactory> DiracEOFactory_;
  RaiiFactoryObj<SolverFactory> SolverFactory_;
  // Objects (Dodwf)
  RaiiFactoryObj<Dirac_DomainWall_EvenOdd> DW5dEO_;
  RaiiFactoryObj<Fopr_DdagD> FoprEO_;
  RaiiFactoryObj<Solver> SolverEO_;
  RaiiFactoryObj<EvenOddUtils::Inverter_WilsonLike> Inv_;
  // Objects (PauliVillars)
  RaiiFactoryObj<Dirac_DomainWall_EvenOdd> DW5dEOpv_;
  RaiiFactoryObj<Fopr_DdagD> FoprEOpv_;
  RaiiFactoryObj<Solver> SolverEOpv_;
  RaiiFactoryObj<EvenOddUtils::Inverter_WilsonLike> InvPV_;
  Dirac_DomainWall_4D* createDirac(InputConfig&);
public:
  DiracDWF4DeoFactory(XML::node node);
};

///////////////////////////////////////////////////////////////////////////////////////////
#ifdef IBM_BGQ_WILSON
/*! @brief Concrete class for creating Dirac_DomainWall_4D_eoSolv with BGQ solv*/
class DiracDWF4dBGQeoFactory : public DiracDWF4dFactory{
  XML::node Dirac_node_;
  // Factories
  RaiiFactoryObj<DiracEvenOdd_DWF5dFactory> DiracEOFactory_;
  RaiiFactoryObj<SolverCG_DWF_opt_Factory> SolverFactory_;
  // Objects (Dodwf)
  RaiiFactoryObj<Dirac_DomainWall_EvenOdd> DW5dEO_;
  RaiiFactoryObj<Solver> SolverEO_;
  RaiiFactoryObj<EvenOddUtils::Inverter_WilsonLike> Inv_;
  // Objects (PauliVillars)
  RaiiFactoryObj<Dirac_DomainWall_EvenOdd> DW5dEOpv_;
  RaiiFactoryObj<Solver> SolverEOpv_;
  RaiiFactoryObj<EvenOddUtils::Inverter_WilsonLike> InvPV_;

  Dirac_DomainWall_4D* createDirac(InputConfig&);
public:
  DiracDWF4dBGQeoFactory(XML::node node);
};
#endif
////////////////////////////////////////////////////////////////////////////////


//////////////
/*! @brief Concrete class for creating Dirac_staggered_EvenOdd */
class DiracStaggeredEvenOddFactory: public DiracStaggeredEvenOddLikeFactory{
  XML::node Dirac_node_;  
  Dirac_staggered_EvenOdd* createDirac(InputConfig&);
  Dirac_staggered_EvenOdd* createDoo(InputConfig&);
public:
  DiracStaggeredEvenOddFactory(const XML::node& node):Dirac_node_(node){}
};

/////////////
/*! @brief Concrete class for creating Dirac_staggered_EvenOdd_Adjoint */
//Valid only in the case of NC=3
#if NC_==3
class DiracStaggeredEvenOddAdjointFactory: public DiracStaggeredEvenOddLikeFactory{
  XML::node Dirac_node_;  
  Dirac_staggered_EvenOdd_Adjoint* createDirac(InputConfig&);
  Dirac_staggered_EvenOdd_Adjoint* createDoo(InputConfig&);
public:
  DiracStaggeredEvenOddAdjointFactory(const XML::node& node):Dirac_node_(node){}
};
#endif

//Add new factories here
//....

#endif
