/*! @file dirac_BFM_wrapper_factory.hpp 
 *  @brief Declaration of BFM Dirac operators factories
 * Time-stamp: <2014-01-30 19:50:14 noaki>
 */
#ifndef DIRAC_BFM_FACT_
#define DIRAC_BFM_FACT_

#include "Dirac_ops/dirac_Operator_Factory.hpp"
#include "Dirac_ops/dirac_Operator_FactoryCreator.hpp"
#include "Solver/solver_Factory.hpp"
#include "dirac_BFM_wrapper.hpp"

//////////////////////////////////
// 5d Operator
//////////////////////////////////
/*! @brief Concrete class for creating Dirac Optimal DWF-5d operators */
class DiracBFMoperatorFactory : public DiracDWF5dFactory {
  const XML::node Dirac_node_;
  // RaiiFactoryObj<Dirac_BFM_Wrapper> BFMop_;

  // Internal operator must be EvenOdd
  RaiiFactoryObj<DiracEvenOdd_DWF5dFactory> Dirac5D_EO_factory_;
  RaiiFactoryObj<DiracWilsonLike_EvenOdd> DiracWilsonEO_;
  Dirac_BFM_Wrapper* createDirac(InputConfig&);
  Dirac_BFM_Wrapper* createDiracPV(InputConfig&);
public:
  DiracBFMoperatorFactory(const XML::node node);
  Dirac_BFM_Wrapper* getDirac(InputConfig& );
  Dirac_BFM_Wrapper* getDiracPV(InputConfig& );
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
class DiracDWF4dBFMeoFactory : public DiracDWF4dFactory{
  const XML::node Dirac_node_;/*< The 4d object level */
  XML::node node_BFM_;/*< The 5d object level */
  XML::node Solver_node_;/*< The Solver object level */
  // Factories
  RaiiFactoryObj<DiracBFMoperatorFactory> DiracBFMFactory_;
  RaiiFactoryObj<SolverCG_DWF_opt_Factory> SolverFactory_;

  // Objects (Dodwf)
  RaiiFactoryObj<Dirac_BFM_Wrapper> BFM_Kernel_;
  RaiiFactoryObj<Solver> SolverEO_;
  RaiiFactoryObj<EvenOddUtils::Inverter_WilsonLike> Inv_;

  // Objects (PauliVillars)
  RaiiFactoryObj<Dirac_BFM_Wrapper> BFM_KernelPV_;
  RaiiFactoryObj<Solver> SolverEOpv_;
  RaiiFactoryObj<EvenOddUtils::Inverter_WilsonLike> InvPV_;

  Dirac_DomainWall_4D* createDirac(InputConfig&);
public:
  DiracDWF4dBFMeoFactory(XML::node node);
  ~DiracDWF4dBFMeoFactory(){
    CCIO::cout << "Destroying DiracDWF4dBFMeoFactory\n";
  };

};

#endif
