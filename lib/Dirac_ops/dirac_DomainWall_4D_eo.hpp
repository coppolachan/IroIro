/*!
 * @file dirac_DomainWall_4D_eo.hpp
 * @brief Definition of Dirac_optimalDomainWall_4D_eo class, whose interface is
 * act like Dirac_optimalDomainWall_4D but contains inverters for Dirac_EvenOdd.
 */
#ifndef DIRAC_OPTIMALDOMAINWALL_4D_EO_INCLUDED
#define DIRAC_OPTIMALDOMAINWALL_4D_EO_INCLUDED

#include "Dirac_ops/dirac_DomainWall.hpp"
#include "Dirac_ops/eoUtils.hpp"
#include "Solver/solver.hpp"
#include "Solver/solver_CG.h"
#include "include/fopr.h"

class Format_F;

class Dirac_optimalDomainWall_4D_eo : public DiracWilsonLike{

  const Dirac_optimalDomainWall D_;            /*!< @brief DWF_5D operator*/
  const Dirac_optimalDomainWall Dpv_;          /*!< @brief DWF_5D operator(PV)*/
  const EvenOddUtils::Iinverter_WilsonLike invD_;  /*!< @brief eo-inverter    */
  const EvenOddUtils::Iinverter_WilsonLike invDpv_;/*!< @brief eo-inverter(PV)*/
  int f5size_;

public:
  Dirac_optimalDomainWall_4D_eo(const Dirac_optimalDomainWall& D,
				const Dirac_optimalDomainWall& Dpv,
				const EvenOddUtils::Inverter_WilsonLike& invD,
				const EvenOddUtils::Inverter_WilsonLike& invDpv)
    :D_(D),Dpv_(Dpv),f5size_(D_.f5size()),
     invD_(invD),invDpv_(invDpv){}

  ~Dirac_optimalDomainWall_4D_eo(){}
  
  size_t fsize() const {return D_.f4size();}
  size_t gsize() const {return D_.gsize(); }

  const Field operator()(int, const Field&) const{}

  const Field mult    (const Field&)const;
  const Field mult_dag(const Field&)const;

  const Field mult_inv    (const Field&)const;
  const Field mult_dag_inv(const Field&)const;

  const Field gamma5  (const Field&f)const{return D_.gamma5_4d(f);}
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
    return Format::Format_F(D_.get_fermionFormat().Nvol());
  }
  double getMass() const {return D_.getMass();}
  const std::vector<int> get_gsite() const{ return D_.get_gsite();}
  const Field md_force(const Field& eta,const Field& zeta) const{}
  void update_internal_state(){}
};

#endif
