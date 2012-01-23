/*!
 * @file dirac_DomainWall_4D.hpp
 *
 * @brief Definition of Dirac_optimalDomainWall_4D class
 *
 */
#ifndef DIRAC_OPTIMALDOMAINWALL_4D_INCLUDED
#define DIRAC_OPTIMALDOMAINWALL_4D_INCLUDED

#include "Dirac_ops/dirac_DomainWall.hpp"
#include "Solver/solver.hpp"
#include "Solver/solver_BiCGStab.h"
#include "include/fopr.h"

class Format_F;

/*!
 * @brief Container for parameter of the 4d Optimal Domain Wall operator
 *
 */
class Dirac_optimalDomainWall_4D : public DiracWilsonLike
{
  const Dirac_optimalDomainWall Dodw_;/*!< @brief 5d Domain Wall fermion operator*/
  const Dirac_optimalDomainWall Dpv_;/*!< @brief 5d Domain Wall fermion operator 
				      *with mass 1.0 (Pauli Villars operator) */
  const Solver* slv_odw_;/*!< @brief %Solver for the Domain Wall fermion operator*/
  const Solver* slv_pv_;/*!< @brief %Solver for the Pauli Villars operator */
public:
  /*!
   * @brief Default constructor using Solver_BiCGStab class as Solver
   */
  Dirac_optimalDomainWall_4D(const Dirac_optimalDomainWall& D,
			     const double stp_cnd_odw, 
			     const double stp_cnd_pv, 
			     const int Niter)
    :Dodw_(D),Dpv_(D,PauliVillars),
     slv_odw_(new Solver_BiCGStab(stp_cnd_odw,Niter,
				  new Fopr_DdagD_Precondition(&Dodw_))),
     slv_pv_(new Solver_BiCGStab(stp_cnd_pv,Niter,
				 new Fopr_DdagD_Precondition(&Dpv_))){}

   /*!
   * @brief Constructor using external solvers (mostly used by factories)
   */
  Dirac_optimalDomainWall_4D(const Dirac_optimalDomainWall& D,
			     const Solver* SolverODWF, 
			     const Solver* SolverPV)
    :Dodw_(D),Dpv_(D,PauliVillars),
     slv_odw_(SolverODWF),slv_pv_(SolverPV){}

  ~Dirac_optimalDomainWall_4D(){}
  
  size_t fsize() const {return Dodw_.f4size();}
  size_t gsize() const {return Dodw_.gsize(); }

  const Field operator()(int, const Field&) const{}

  const Field mult    (const Field&)const;
  const Field mult_dag(const Field&)const;
  const Field gamma5  (const Field&f)const{return Dodw_.gamma5_4d(f);}

  const Field mult_inv    (const Field&)const;
  const Field mult_dag_inv(const Field&)const;

  ////////////////////////////////////////Preconditioned versions
  // 4d operator has no preconditioner now 
  const Field mult_prec     (const Field&f)const{return f;}
  const Field mult_dag_prec (const Field&f)const{return f;}
  const Field left_prec     (const Field&f)const{return f;}
  const Field right_prec    (const Field&f)const{return f;}
  const Field left_dag_prec (const Field&f)const{return f;}
  const Field right_dag_prec(const Field&f)const{return f;}
  //////////////////////////////////////////////////////////////

  const Field signKernel(const Field&)const;

  const Format::Format_F get_fermionFormat() const {
    return Dodw_.get_fermionFormat();
  }
  const double getMass() const {return Dodw_.getMass();}
  const std::vector<int> get_gsite() const{ return Dodw_.get_gsite();}
  const Field md_force(const Field& eta,const Field& zeta) const{}
  
  void update_internal_state(){}
};

#endif
