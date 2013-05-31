/*!
 * @file dirac_DomainWall_4D_fullSolv.hpp
 * @brief Definition of Dirac_optimalDomainWall_4D class with full Solver
 Time-stamp: <2013-05-31 07:16:04 noaki>
 */
#ifndef DIRAC_OPTIMALDOMAINWALL_4D_FULLSOLV_INCLUDED
#define DIRAC_OPTIMALDOMAINWALL_4D_FULLSOLV_INCLUDED

#include "Dirac_ops/dirac_DomainWall.hpp"
#include "Solver/solver.hpp"
#include "Solver/solver_CG.hpp"
#include "include/fopr.h"
#include "EigenModes/lowModesHandler.hpp"
/*!
 * @brief Container for parameter of the 4d Optimal Domain Wall operator
 */
class Dirac_optimalDomainWall_4D_fullSolv : public Dirac_optimalDomainWall_4D{
  const Dirac_optimalDomainWall* Dodw_;/*!< @brief 5d Domain Wall fermion operator*/
  const Dirac_optimalDomainWall* Dpv_;/*!< @brief 5d Domain Wall fermion operator 
				      *with mass 1.0 (Pauli Villars operator) */
  const Solver* slv_odw_;/*!< @brief %Solver for the Domain Wall fermion operator*/
  const Solver* slv_pv_;/*!< @brief %Solver for the Pauli Villars operator */
  double mq_;
  const LowModesHandler* lmh_;    /*!< @brief LowModes of Kernel operator */ 
  bool lmhGiven_;
  size_t fsize_;
  void mult_normal(Field&,const Field&)const;
  void mult_lmp(Field&,const Field&)const;
  void mult_inv_normal(Field&,const Field&)const;
  void mult_inv_lmp(Field&,const Field&)const;
  void(Dirac_optimalDomainWall_4D_fullSolv::*mult_core)(Field&,const Field&)const;
  void(Dirac_optimalDomainWall_4D_fullSolv::*mult_inv_core)(Field&,const Field&)const;
public:
   /*!
   * @brief Constructor using external solvers (mostly used by factories)
   */
  Dirac_optimalDomainWall_4D_fullSolv(const Dirac_optimalDomainWall* D,
				      const Dirac_optimalDomainWall* Dpv,
				      const Solver* SolverODWF, 
				      const Solver* SolverPV,
				      const LowModesHandler* lmh = NULL)
    :Dodw_(D),Dpv_(Dpv),slv_odw_(SolverODWF),slv_pv_(SolverPV),
     mq_(D->getMass()),lmh_(lmh),lmhGiven_(false),fsize_(D->fsize()),
     mult_core(&Dirac_optimalDomainWall_4D_fullSolv::mult_normal),
     mult_inv_core(&Dirac_optimalDomainWall_4D_fullSolv::mult_inv_normal){
    if(lmh_){
      lmhGiven_= true;
      mult_core = &Dirac_optimalDomainWall_4D_fullSolv::mult_lmp;
      mult_inv_core = &Dirac_optimalDomainWall_4D_fullSolv::mult_inv_lmp;
    }
  }

  Dirac_optimalDomainWall_4D_fullSolv(XML::node node,
				      const Dirac_optimalDomainWall* D,
				      const Dirac_optimalDomainWall* Dpv,
				      const Solver* SolverODWF, 
				      const Solver* SolverPV)
    :Dodw_(D),Dpv_(Dpv),slv_odw_(SolverODWF),slv_pv_(SolverPV),
     mq_(D->getMass()),lmh_(NULL),lmhGiven_(false),fsize_(D->fsize()),
     mult_core(&Dirac_optimalDomainWall_4D_fullSolv::mult_normal),
     mult_inv_core(&Dirac_optimalDomainWall_4D_fullSolv::mult_inv_normal){

    XML::node LMPnode = node.child("LowModesPrecondition");
    if(LMPnode !=NULL){
      lmh_= new LowModesHandler(LMPnode);
      mult_core = &Dirac_optimalDomainWall_4D_fullSolv::mult_lmp;
      mult_inv_core = &Dirac_optimalDomainWall_4D_fullSolv::mult_inv_lmp;
    }
    XML::descend(node,"Kernel5d",MANDATORY);
  }

  ~Dirac_optimalDomainWall_4D_fullSolv(){ if(lmh_ && !lmhGiven_) delete lmh_;}
  
  size_t fsize() const {return Dodw_->f4size();}
  size_t gsize() const {return Dodw_->gsize(); }

  const Field mult    (const Field&)const;
  const Field mult_dag(const Field&)const;
  const Field gamma5  (const Field&)const;

  const Field mult_inv    (const Field&)const;
  const Field mult_dag_inv(const Field&)const;

  const Field signKernel(const Field&)const;
  double getMass() const {return mq_;}
};

#endif
