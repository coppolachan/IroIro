/*!
 * @file dirac_DomainWall_4D_fullSolv.hpp
 * @brief Definition of Dirac_DomainWall_4D class with full Solver
 Time-stamp: <2014-07-02 21:31:57 noaki>
 */
#ifndef DIRAC_DOMAINWALL_4D_FULLSOLV_INCLUDED
#define DIRAC_DOMAINWALL_4D_FULLSOLV_INCLUDED

#include "dirac_DomainWall.hpp"
#include "Solver/solver.hpp"
#include "Solver/solver_CG.hpp"
#include "Fopr/fopr.h"

enum DW5dPrecond{NoPrecond,LUprecond};

/*!
 * @brief Container for parameter of the 4d Optimal Domain Wall operator
 */
class Dirac_DomainWall_4D_fullSolv : public Dirac_DomainWall_4D{
  const Dirac_DomainWall* Dodw_;/*!< @brief 5d Domain Wall fermion operator*/
  const Dirac_DomainWall* Dpv_;/*!< @brief 5d Domain Wall fermion operator 
				      *with mass 1.0 (Pauli Villars operator) */
  const Solver* slv_odw_;/*!< @brief %Solver for the Domain Wall fermion op*/
  const Solver* slv_pv_;/*!< @brief %Solver for the Pauli Villars operator */

  const Field Bproj(const Field&) const;
  const Field Bproj_dag(const Field&) const;


  void(Dirac_DomainWall_4D_fullSolv::*mult_core)(Field&,const Field&)const;
  void(Dirac_DomainWall_4D_fullSolv::*mult_inv_core)(Field&,const Field&)const;

  void mult_std(    Field&,const Field&)const;
  void mult_inv_std(Field&,const Field&)const;

  void mult_LU(    Field&,const Field&)const;
  void mult_inv_LU(Field&,const Field&)const;

  DW5dMatrix d5_;
  double mq_;
  size_t fsize_;
public:
   /*!
   * @brief Constructor using external solvers (mostly used by factories)
   */
  Dirac_DomainWall_4D_fullSolv(const Dirac_DomainWall* D,
			       const Dirac_DomainWall* Dpv,
			       const Solver* SolverDW, 
			       const Solver* SolverPV,
			       DW5dPrecond precond=NoPrecond)
    :Dodw_(D),Dpv_(Dpv),slv_odw_(SolverDW),slv_pv_(SolverPV),
     mult_core(    &Dirac_DomainWall_4D_fullSolv::mult_std),
     mult_inv_core(&Dirac_DomainWall_4D_fullSolv::mult_inv_std),
     mq_(D->getMass()),fsize_(D->f4size()){
    if(precond==LUprecond){
      mult_core =     &Dirac_DomainWall_4D_fullSolv::mult_LU;
      mult_inv_core = &Dirac_DomainWall_4D_fullSolv::mult_inv_LU;
    }
  }
  
  size_t fsize() const {return fsize_;}
  size_t gsize() const {return Dodw_->gsize(); }
  double getMass() const {return mq_;}

  const Field* getGaugeField_ptr()const{ return Dodw_->getGaugeField_ptr(); }

  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;
  const Field mult_inv(const Field&)const;
  const Field mult_dag_inv(const Field&)const;

  const Field gamma5(const Field&)const;
};

#endif
