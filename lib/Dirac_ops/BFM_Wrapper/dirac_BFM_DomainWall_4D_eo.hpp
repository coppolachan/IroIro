/*!
 * @file dirac_BFM_DomainWall_4D_eo.hpp
 * @brief Definition of Dirac_BFM_DomainWall_4D_eo class
 *
 * Time-stamp: <2014-07-10 18:41:22 noaki>
 */
#ifndef DIRAC_BFM_DWF_4D_EO_INCLUDED
#define DIRAC_BFM_DWF_4D_EO_INCLUDED

#include "Dirac_ops/dirac_DomainWall.hpp"
#include "Dirac_ops/eoUtils.hpp"
#include "dirac_BFM_wrapper.hpp"
#include "Solver/solver.hpp"
#include "Solver/solver_CG.hpp"
#include "Fopr/fopr.h"

class Dirac_BFM_DomainWall_4D_eo : public Dirac_DomainWall_4D{
private:
  const EvenOddUtils::Inverter_WilsonLike* invD_;  /*!< @brief eo-inverter    */
  const EvenOddUtils::Inverter_WilsonLike* invDpv_;/*!< @brief eo-inverter(PV)*/

  const Field Bproj(const Field&)const;
  const Field Bproj_dag(const Field&)const;

  Dirac_BFM_Wrapper* BFM_Op_;/*!< @brief 5D DWF operator with mass m */
  Dirac_BFM_Wrapper* BFM_Op_PV_;/*!< @brief 5D DWF operator with mass 1, Pauli-Villars */


  int Nvol_,N5_;
  ffmt_t ff_;
  size_t fsize_;  
  double mq_;
  GammaMatrix dm_;
  DW5dMatrix d5_;
public:
  Dirac_BFM_DomainWall_4D_eo(XML::node node,
			     const EvenOddUtils::Inverter_WilsonLike* invD,
			     const EvenOddUtils::Inverter_WilsonLike* invDpv,
			     Dirac_BFM_Wrapper* BFM_DWF,
			     Dirac_BFM_Wrapper* BFM_DWFpv)
    :BFM_Op_(BFM_DWF),BFM_Op_PV_(BFM_DWFpv),
     invD_(invD),invDpv_(invDpv),
     Nvol_(CommonPrms::instance()->Nvol()),
     ff_(Nvol_),fsize_(ff_.size()){
    XML::read(node,"N5d",N5_,MANDATORY);
    XML::read(node,"mass",mq_,MANDATORY);
  }
  
  Dirac_BFM_DomainWall_4D_eo(int N5,double mq,
			     const EvenOddUtils::Inverter_WilsonLike* invD,
			     const EvenOddUtils::Inverter_WilsonLike* invDpv)
    :invD_(invD),invDpv_(invDpv),
     Nvol_(CommonPrms::instance()->Nvol()),
     ff_(Nvol_),fsize_(ff_.size()),N5_(N5),mq_(mq){}

  ~Dirac_BFM_DomainWall_4D_eo();

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
