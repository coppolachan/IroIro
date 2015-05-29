/*! @file dirac_BFM_wrapper_factory.hpp 
 *  @brief Declaration of BFM Dirac operators factories
 *
 * Time-stamp: <2015-05-18 15:47:18 cossu>
 */
#ifndef DIRAC_BFM_FACT_
#define DIRAC_BFM_FACT_

#include "Dirac_ops/dirac_Operator_Factory.hpp"
#include "Dirac_ops/dirac_Operator_FactoryCreator.hpp"
#include "Solver/solver_Factory.hpp"
#include "dirac_BFM_wrapper.hpp"
#include "dirac_BFM_DomainWall_4D_eo.hpp"

//////////////////////////////////
// 5d Operator
//////////////////////////////////
/*! @brief Concrete class for creating Dirac DWF-5d operator with BFM */
class DiracBFMoperatorFactory : public DiracDWF5dFactory {
  const XML::node Dirac_node_;
 
  // Internal operator must be EvenOdd
  RaiiFactoryObj<DiracEvenOdd_DWF5dFactory> Dirac5D_EO_factory_;
  RaiiFactoryObj<DiracWilsonLike_EvenOdd> DiracWilsonEO_;
  Dirac_BFM_Wrapper* createDirac(const InputConfig&);
  Dirac_BFM_Wrapper* createDiracPV(const InputConfig&);
public:
  DiracBFMoperatorFactory(const XML::node node);
  Dirac_BFM_Wrapper* getDirac(const InputConfig& );
  Dirac_BFM_Wrapper* getDiracPV(const InputConfig& );
};


//////////////////////////////
// 4d Operator
//////////////////////////////

/*! @brief Concrete class for creating Dirac_optimalDomainWall_4D_eoSolv with BFM solv

  XML structure
  operator name: DiracOptimalDomainWall4d_BFMeo

  internal structure (similar to the action construction for BFM
  < .... name="DiracOptimalDomainWall4d_BFMeo">
    <Kernel5d>     comment: no name needed since it creates a BFM Object
       <Operator>   comment: no name needed, it is a 5D domain wall operator E/O
          definition of the 5D DWF operator with usual rules
       </Operator>

       <BFMKernel>
         <mass>...</mass>
	 <M5>...</M5>
	 <N5d>...</N5d>
	 <scale>...</scale>
       </BFMKernel>
    <Kernel5d>

    <Solver_DWF-EO_BGQ>
        usual definition of the solver
    </Solver_DWF-EO_BGQ>
 */

/* templated class to allow to the HDCG operator */
template < class BFM5dFact >
class DiracDWF4dBFMeoFactory  : public DiracDWF4dFactory{
  const XML::node Dirac_node_;/*< The 4d object level */
  XML::node node_BFM_;/*< The 5d object level */
  XML::node Solver_node_;/*< The Solver object level */
  // Factories
  RaiiFactoryObj< BFM5dFact > DiracBFMFactory_;
  RaiiFactoryObj< DiracBFMoperatorFactory > DiracBFMFactory_PV_; // differentiate PV in the general case
  RaiiFactoryObj<SolverCG_DWF_opt_Factory> SolverFactory_;

  // Objects (Dodwf)
  RaiiFactoryObj<Dirac_BFM_Wrapper> BFM_Kernel_;
  RaiiFactoryObj<Solver> SolverEO_;
  RaiiFactoryObj<EvenOddUtils::Inverter_WilsonLike> Inv_;

  // Objects (PauliVillars)
  RaiiFactoryObj<Dirac_BFM_Wrapper> BFM_KernelPV_;
  RaiiFactoryObj<Solver> SolverEOpv_;
  RaiiFactoryObj<EvenOddUtils::Inverter_WilsonLike> InvPV_;

  Dirac_DomainWall_4D* createDirac(const InputConfig&);
public:
  DiracDWF4dBFMeoFactory(XML::node node);
  ~DiracDWF4dBFMeoFactory(){};

};


///////////////////// 4D Operator templated methods definition
template < class BFM5dFact > 
DiracDWF4dBFMeoFactory<BFM5dFact>::DiracDWF4dBFMeoFactory(XML::node node)
:Dirac_node_(node){
  node_BFM_ = Dirac_node_;
  XML::descend(node_BFM_,"Kernel5d",MANDATORY);
  DiracBFMFactory_.save(new BFM5dFact(node_BFM_));
  DiracBFMFactory_PV_.save(new DiracBFMoperatorFactory(node_BFM_));



  Solver_node_ = node;
  XML::descend(Solver_node_,"Solver_DWF-EO_BGQ",MANDATORY);
  SolverFactory_.save(new SolverCG_DWF_opt_Factory(Solver_node_));
}

template < class BFM5dFact > 
Dirac_DomainWall_4D* DiracDWF4dBFMeoFactory<BFM5dFact>::createDirac(const InputConfig& input){
  BFM_KernelPV_.save( DiracBFMFactory_PV_.get()->getDiracPV(input)); 
  BFM_Kernel_.save( DiracBFMFactory_.get()->getDirac(input));
 

  SolverEO_.save(SolverFactory_.get()->getSolver(BFM_Kernel_.get()));
  SolverEOpv_.save(SolverFactory_.get()->getSolver(BFM_KernelPV_.get()));

  BFM_Kernel_.get()->set_SolverParams(Solver_node_);
  BFM_KernelPV_.get()->set_SolverParams(Solver_node_);
  
  BFM_Kernel_.get()->initialize(node_BFM_);
  BFM_KernelPV_.get()->initialize();

  Inv_.save(new EvenOddUtils::Inverter_WilsonLike(BFM_Kernel_.get()->getInternalEO(),SolverEO_.get()));
  InvPV_.save(new EvenOddUtils::Inverter_WilsonLike(BFM_KernelPV_.get()->getInternalEO(),SolverEOpv_.get()));

  XML::node current_node = node_BFM_;  
  XML::descend(current_node, "Operator", MANDATORY);  
  return new Dirac_BFM_DomainWall_4D_eo(current_node,
					Inv_.get(), 
					InvPV_.get(),
					BFM_Kernel_.get(),
					BFM_KernelPV_.get());
}



#endif
