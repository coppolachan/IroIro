/*!
 * @file dirac_DomainWall_4D_eoSolv.hpp
 * @brief Definition of Dirac_optimalDomainWall_4D_eoSolv class
 * which contains e/o-solver
 Time-stamp: <2013-05-02 17:23:09 noaki>
 */
#ifndef DIRAC_OPTIMALDOMAINWALL_4D_EOSOLV_INCLUDED
#define DIRAC_OPTIMALDOMAINWALL_4D_EOSOLV_INCLUDED

#include "Dirac_ops/dirac_DomainWall.hpp"
#include "Dirac_ops/eoUtils.hpp"
#include "Solver/solver.hpp"
#include "Solver/solver_CG.hpp"
#include "include/fopr.h"

class Dirac_optimalDomainWall_4D_eoSolv : public Dirac_optimalDomainWall_4D{

  const Dirac_optimalDomainWall* D_;        /*!< @brief DWF_5D operator*/
  const Dirac_optimalDomainWall* Dpv_;      /*!< @brief DWF_5D operator(PV)*/
  const EvenOddUtils::Inverter_WilsonLike* invD_;  /*!< @brief eo-inverter    */
  const EvenOddUtils::Inverter_WilsonLike* invDpv_;/*!< @brief eo-inverter(PV)*/

public:
  Dirac_optimalDomainWall_4D_eoSolv(const Dirac_optimalDomainWall* D,
				    const Dirac_optimalDomainWall* Dpv,
				    const EvenOddUtils::Inverter_WilsonLike* invD,
				    const EvenOddUtils::Inverter_WilsonLike* invDpv)
    :D_(D),Dpv_(Dpv),invD_(invD),invDpv_(invDpv){}

  ~Dirac_optimalDomainWall_4D_eoSolv(){}
  
  size_t fsize() const {return D_->f4size();}
  size_t gsize() const {return D_->gsize(); }

  const Field mult    (const Field&)const;
  const Field mult_dag(const Field&)const;
  const Field gamma5  (const Field&)const;

  const Field mult_inv    (const Field&)const;
  const Field mult_dag_inv(const Field&)const;

  const Field signKernel(const Field&)const;
  double getMass() const {return D_->getMass();}
};

#endif
