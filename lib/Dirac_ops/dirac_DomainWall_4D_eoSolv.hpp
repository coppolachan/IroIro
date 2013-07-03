/*!
 * @file dirac_DomainWall_4D_eoSolv.hpp
 * @brief Definition of Dirac_optimalDomainWall_4D_eoSolv class
 * which contains e/o-solver
 Time-stamp: <2013-07-04 02:16:34 noaki>
 */
#ifndef DIRAC_OPTIMALDOMAINWALL_4D_EOSOLV_INCLUDED
#define DIRAC_OPTIMALDOMAINWALL_4D_EOSOLV_INCLUDED

#include "dirac_DomainWall.hpp"
#include "eoUtils.hpp"
#include "Solver/solver.hpp"
#include "Solver/solver_CG.hpp"
#include "include/fopr.h"

class Dirac_optimalDomainWall_4D_eoSolv : public Dirac_optimalDomainWall_4D{
private:
  const EvenOddUtils::Inverter_WilsonLike* invD_;  /*!< @brief eo-inverter    */
  const EvenOddUtils::Inverter_WilsonLike* invDpv_;/*!< @brief eo-inverter(PV)*/

  const Field Bproj(const Field&)const;
  const Field Bproj_dag(const Field&)const;

  int Nvol_,N5_;
  ffmt_t ff_;
  size_t fsize_;  
  double mq_;

  // For ExactLowModePrecond
public:
  Dirac_optimalDomainWall_4D_eoSolv(XML::node node,
				    const EvenOddUtils::Inverter_WilsonLike* invD,
				    const EvenOddUtils::Inverter_WilsonLike* invDpv)
    :invD_(invD),invDpv_(invDpv),
     Nvol_(CommonPrms::instance()->Nvol()),
     ff_(Nvol_),fsize_(ff_.size()){
    XML::read(node,"N5d",N5_,MANDATORY);
    XML::read(node,"mass",mq_,MANDATORY);
  }
  
  Dirac_optimalDomainWall_4D_eoSolv(int N5,double mq,
				    const EvenOddUtils::Inverter_WilsonLike* invD,
				    const EvenOddUtils::Inverter_WilsonLike* invDpv)
    :invD_(invD),invDpv_(invDpv),
     Nvol_(CommonPrms::instance()->Nvol()),N5_(N5),ff_(Nvol_),
     fsize_(ff_.size()),mq_(mq){}

  size_t fsize() const {return fsize_;}
  size_t gsize()const{
    return Nvol_*gfmt_t::Nin()*CommonPrms::instance()->Ndim(); }
  double getMass() const {return mq_;}

  const Field* getGaugeField_ptr()const{ return invD_->getGaugeField_ptr(); }

  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;
  const Field mult_inv(const Field&)const;
  const Field mult_dag_inv(const Field&)const;

  const Field gamma5  (const Field&)const;
};

#endif
