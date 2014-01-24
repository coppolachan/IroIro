/*!
 * @file dirac_DomainWall_adjoint.hpp
 * @brief Declaration of class Dirac_DomainWall_Adjoint (5d operator)
 Time-stamp: <2014-01-24 14:44:51 neo>
 */
#ifndef DIRAC_DOMAINWALL_ADJOINT_INCLUDED
#define DIRAC_DOMAINWALL_ADJOINT_INCLUDED

#include "domainWallCore.hpp"

class Dirac_DomainWall_Adjoint : public DiracWilsonLike {
private:
  DomainWallParams prms_;
  ProjPadj projP_;
  ProjMadj projM_;
  Op5Dadj op5D_;

  const DiracWilsonLike* Dw_;
  DomainWallCore DWcore_;
public:
  ////////////// Constructors for full indexing //////////////////
  /*! @brief xml constructor */
  Dirac_DomainWall_Adjoint(XML::node dw_node,const DiracWilsonLike* Dw,
			   DWFType Type= Regular)
    :prms_(dw_node,Type),Dw_(Dw),
     projP_(CommonPrms::instance()->Nvol(),prms_.N5_),
     projM_(CommonPrms::instance()->Nvol(),prms_.N5_),
     op5D_( CommonPrms::instance()->Nvol(),prms_.N5_),
     DWcore_(prms_,Dw,projP_,projM_){
    //
    assert(Dw->fsize() == afmt_t::Nin()*CommonPrms::instance()->Nvol());}
  /*! @brief  manual constructor */
  Dirac_DomainWall_Adjoint(double b,double c,double M0,double mq,
			   const std::vector<double>& omega,
			   const DiracWilsonLike* Dw)
    :prms_(b,c,M0,mq,omega),Dw_(Dw),
     projP_(CommonPrms::instance()->Nvol(),prms_.N5_),
     projM_(CommonPrms::instance()->Nvol(),prms_.N5_),
     op5D_(CommonPrms::instance()->Nvol(),prms_.N5_),
     DWcore_(prms_,Dw,projP_,projM_){
    //
    assert(Dw->fsize() == afmt_t::Nin()*CommonPrms::instance()->Nvol());}

  /*! @brief copy constructor*/
  Dirac_DomainWall_Adjoint(const Dirac_DomainWall_Adjoint& Dc, 
			   DWFType Type=Regular)
    :prms_(Dc.prms_,Type),Dw_(Dc.Dw_),
     projP_(Dc.projP_),projM_(Dc.projM_),op5D_(Dc.op5D_),
     DWcore_(prms_,Dw_,projP_,projM_){}

  ////////////// Constructors for e/o indexing //////////////////
  /*! @brief xml constructor */
  Dirac_DomainWall_Adjoint(XML::node dw_node,const DiracWilsonLike* Dw,
			   DWFutils::EvenOdd_tag EO,DWFType Type= Regular)
    :prms_(dw_node,Type),Dw_(Dw),
     projP_(CommonPrms::instance()->Nvol()/2,prms_.N5_),
     projM_(CommonPrms::instance()->Nvol()/2,prms_.N5_),
     op5D_(CommonPrms::instance()->Nvol()/2,prms_.N5_),
     DWcore_(prms_,Dw,projP_,projM_,EO){
    //
    assert(Dw->fsize() == afmt_t::Nin()*CommonPrms::instance()->Nvol()/2);}

  /*! @brief manual constructor*/
  Dirac_DomainWall_Adjoint(double b,double c,double M0,double mq,
			   const std::vector<double>& omega,
			   const DiracWilsonLike* Dw,
			   DWFutils::EvenOdd_tag EO)
    :prms_(b,c,M0,mq,omega),Dw_(Dw),
     projP_(CommonPrms::instance()->Nvol()/2,prms_.N5_),
     projM_(CommonPrms::instance()->Nvol()/2,prms_.N5_),
     op5D_(CommonPrms::instance()->Nvol()/2,prms_.N5_),
     DWcore_(prms_,Dw,projP_,projM_,EO){
    //
    assert(Dw->fsize() == afmt_t::Nin()*CommonPrms::instance()->Nvol()/2);}

  /*! @brief copy constructor*/
  Dirac_DomainWall_Adjoint(const Dirac_DomainWall_Adjoint& Dc, 
			   DWFutils::EvenOdd_tag EO,DWFType Type=Regular)
    :prms_(Dc.prms_,Type),Dw_(Dc.Dw_),
     projP_(Dc.projP_),projM_(Dc.projM_),op5D_(Dc.op5D_),
     DWcore_(prms_,Dw_,projP_,projM_,EO){}
  ////////////////

  size_t fsize() const{ return Dw_->fsize()*prms_.N5_;}
  size_t f4size() const{ return Dw_->fsize();}
  size_t gsize()const{return Dw_->gsize();}
  double getMass() const{return prms_.mq_;}
  int getN5() const{return prms_.N5_;}

  const Field* getGaugeField_ptr()const{ return Dw_->getGaugeField_ptr(); }

  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;
  const Field gamma5(const Field&)const;

  /*@! brief mult in the heavy M0 limit */
  const Field mult_hop5(const Field& f5)const;    
  const Field mult_hop5_inv(const Field& f5)const;
  const Field mult_hop5_dag(const Field& f5)const;
  const Field mult_hop5_dinv(const Field& f5)const;

  void gamma5_4d(Field&,const Field&)const;
  void md_force_p(Field&,const Field&,const Field&)const;
  void md_force_m(Field&,const Field&,const Field&)const;

  const Field md_force(const Field& eta,const Field& zeta)const;
  void update_internal_state(){}
};

inline const Field Dirac_DomainWall_Adjoint::
mult(const Field& f5) const{ return DWcore_.mult(f5);}

inline const Field Dirac_DomainWall_Adjoint::
mult_dag(const Field& f5) const{ return DWcore_.mult_dag(f5);}

inline const Field Dirac_DomainWall_Adjoint::
gamma5(const Field& f5) const{ return op5D_.gamma5(f5);}

inline const Field Dirac_DomainWall_Adjoint::
md_force(const Field& eta,const Field& zeta)const{
  return DWcore_.md_force(eta,zeta);}

inline void Dirac_DomainWall_Adjoint::
gamma5_4d(Field& w,const Field& f)const{ w = Dw_->gamma5(f);}

inline const Field Dirac_DomainWall_Adjoint::
mult_hop5(const Field& f5)const{ return DWcore_.mult_hop5(f5);}

inline const Field Dirac_DomainWall_Adjoint::
mult_hop5_dag(const Field& f5)const{ return DWcore_.mult_hop5_dag(f5);}

inline const Field Dirac_DomainWall_Adjoint::
mult_hop5_inv(const Field& f5)const{ return DWcore_.mult_hop5_inv(f5);}

inline const Field Dirac_DomainWall_Adjoint::
mult_hop5_dinv(const Field& f5)const{ return DWcore_.mult_hop5_dinv(f5);}

inline void Dirac_DomainWall_Adjoint::
md_force_p(Field& fce,const Field& phi,const Field& psi)const{
  DWcore_.md_force_p(fce,phi,psi);}  

inline void Dirac_DomainWall_Adjoint::
md_force_m(Field& fce,const Field& phi,const Field& psi)const{
  DWcore_.md_force_m(fce,phi,psi);}  

#endif
