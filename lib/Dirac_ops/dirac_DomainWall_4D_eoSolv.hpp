/*!
 * @file dirac_DomainWall_4D_eoSolv.hpp
 * @brief Definition of Dirac_optimalDomainWall_4D_eoSolv class
 * which contains e/o-solver
 Time-stamp: <2013-05-30 20:59:09 noaki>
 */
#ifndef DIRAC_OPTIMALDOMAINWALL_4D_EOSOLV_INCLUDED
#define DIRAC_OPTIMALDOMAINWALL_4D_EOSOLV_INCLUDED

#include "Dirac_ops/dirac_DomainWall.hpp"
#include "Dirac_ops/eoUtils.hpp"
#include "Solver/solver.hpp"
#include "Solver/solver_CG.hpp"
#include "include/fopr.h"
#include "EigenModes/lowModesHandler.hpp"

class Dirac_optimalDomainWall_4D_eoSolv : public Dirac_optimalDomainWall_4D{
private:
  const EvenOddUtils::Inverter_WilsonLike* invD_;  /*!< @brief eo-inverter    */
  const EvenOddUtils::Inverter_WilsonLike* invDpv_;/*!< @brief eo-inverter(PV)*/

  int Nvol_,N5_;
  ffmt_t ff_;
  size_t fsize_;  
  double mq_;

  // For ExactLowModePrecond
  const LowModesHandler* lmh_;    /*!< @brief LowModes of Kernel operator */ 
  bool lmhGiven_;
  void mult_normal(Field&,const Field&)const;
  void mult_lmp(Field&,const Field&)const;
  void mult_inv_normal(Field&,const Field&)const;
  void mult_inv_lmp(Field&,const Field&)const;
  void(Dirac_optimalDomainWall_4D_eoSolv::*mult_core)(Field&,const Field&)const;
  void(Dirac_optimalDomainWall_4D_eoSolv::*mult_inv_core)(Field&,const Field&)const;

public:
  Dirac_optimalDomainWall_4D_eoSolv(XML::node node,
				    const EvenOddUtils::Inverter_WilsonLike* invD,
				    const EvenOddUtils::Inverter_WilsonLike* invDpv)
    :invD_(invD),invDpv_(invDpv),
     Nvol_(CommonPrms::instance()->Nvol()),
     ff_(Nvol_),fsize_(ff_.size()),lmh_(NULL),lmhGiven_(false),
     mult_core(&Dirac_optimalDomainWall_4D_eoSolv::mult_normal),
     mult_inv_core(&Dirac_optimalDomainWall_4D_eoSolv::mult_inv_normal){

    XML::node LMPnode = node.child("LowModesPrecondition");
    if(LMPnode !=NULL){
      lmh_= new LowModesHandler(LMPnode);
      mult_core = &Dirac_optimalDomainWall_4D_eoSolv::mult_lmp;
      mult_inv_core = &Dirac_optimalDomainWall_4D_eoSolv::mult_inv_lmp;
    }
    XML::descend(node,"Kernel5d",MANDATORY);
    XML::read(node,"N5d",N5_,MANDATORY);
    XML::read(node,"mass",mq_,MANDATORY);
  }

  Dirac_optimalDomainWall_4D_eoSolv(int N5,double mq,
				    const EvenOddUtils::Inverter_WilsonLike* invD,
				    const EvenOddUtils::Inverter_WilsonLike* invDpv,
				    const LowModesHandler* lmh = NULL)
    :invD_(invD),invDpv_(invDpv),
     Nvol_(CommonPrms::instance()->Nvol()),N5_(N5),ff_(Nvol_),
     fsize_(ff_.size()),mq_(mq),lmh_(lmh),lmhGiven_(false),
     mult_core(&Dirac_optimalDomainWall_4D_eoSolv::mult_normal),
     mult_inv_core(&Dirac_optimalDomainWall_4D_eoSolv::mult_inv_normal){
    
    if(lmh_){
      lmhGiven_= true;
      mult_core = &Dirac_optimalDomainWall_4D_eoSolv::mult_lmp;
      mult_inv_core = &Dirac_optimalDomainWall_4D_eoSolv::mult_inv_lmp;
    }
  }  

  ~Dirac_optimalDomainWall_4D_eoSolv(){ if(lmh_ && !lmhGiven_) delete lmh_;}

  size_t fsize() const {return fsize_;}
  size_t gsize()const{
    return Nvol_*gfmt_t::Nin()*CommonPrms::instance()->Ndim(); }

  const Field mult    (const Field&)const;
  const Field mult_dag(const Field&)const;
  const Field gamma5  (const Field&)const;
  const Field Bproj(const Field&)const;
  const Field Bproj_dag(const Field&)const;

  const Field mult_inv    (const Field&)const;
  const Field mult_dag_inv(const Field&)const;

  const Field signKernel(const Field&)const;
  double getMass() const {return mq_;}
};

#endif
