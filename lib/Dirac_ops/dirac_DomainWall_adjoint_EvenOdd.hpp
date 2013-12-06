/*!
 * @file dirac_DomainWall_adjoint_EvenOdd.hpp
 * @brief Declaration of class Dirac_DomainWall_Adjoint_EvenOdd (5d operator)
 Time-stamp: <2013-12-05 11:08:41 noaki>
 */
#ifndef DIRAC_DOMAINWALL_ADJOINT_EVENODD_INCLUDED
#define DIRAC_DOMAINWALL_ADJOINT_EVENODD_INCLUDED

#include "dirac_DomainWall_adjoint.hpp"
/*!
 * @brief Defines the 5d DW-adjoint operator with even/odd site indexing
 */
class Dirac_DomainWall_Adjoint_EvenOdd :public DiracWilsonLike_EvenOdd {
private:
  const Dirac_DomainWall_Adjoint Deo_;
  const Dirac_DomainWall_Adjoint Doe_;

  void md_force_eo(Field&,const Field&,const Field&)const;
  void md_force_oe(Field&,const Field&,const Field&)const;
  
  Dirac_DomainWall_Adjoint_EvenOdd(const Dirac_DomainWall_Adjoint_EvenOdd&);
  /*!< simple copy is prohibited */

public:
  Dirac_DomainWall_Adjoint_EvenOdd(XML::node dw_node,
				   DiracWilsonLike_EvenOdd* Kernel,
				   DWFType Type=Regular)
    :Deo_(dw_node,Kernel->getDeo(),DWF::EvenOdd_tag(),Type),
     Doe_(dw_node,Kernel->getDoe(),DWF::EvenOdd_tag(),Type){
    //
    assert(Kernel->fsize() == afmt_t::Nin()*CommonPrms::instance()->Nvol()/2);
#if VERBOSITY>DEBUG_VERB_LEVEL 
    CCIO::cout<<"Dirac_DomainWall_Evenodd created"<<std::endl;
#endif
  }

  /*! @brief copy constractor to create Pauli-Villars operator */
  Dirac_DomainWall_Adjoint_EvenOdd(const Dirac_DomainWall_Adjoint_EvenOdd& D, 
				   DWFType Type=Regular)
    :Deo_(D.Deo_,DWF::EvenOdd_tag(),Type),
     Doe_(D.Doe_,DWF::EvenOdd_tag(),Type){}
  
  Dirac_DomainWall_Adjoint_EvenOdd(double b,double c,double M0,double mq,
				   const std::vector<double>& omega,
				   const DiracWilsonLike* Keo,
				   const DiracWilsonLike* Koe)
    :Deo_(b,c,M0,mq,omega,Keo,DWF::EvenOdd_tag()),
     Doe_(b,c,M0,mq,omega,Koe,DWF::EvenOdd_tag()){}
  
  ~Dirac_DomainWall_Adjoint_EvenOdd(){}
  
  size_t f4size()const{ return Deo_.f4size();}
  size_t fsize()const{ return Deo_.fsize();}
  size_t gsize()const{ return Deo_.gsize();}

  const DiracWilsonLike* getDeo()const{ return &Deo_;}
  const DiracWilsonLike* getDoe()const{ return &Doe_;}
  
  double getMass() const{return Deo_.getMass();}
  int getN5() const{return Deo_.getN5();}

  const Field* getGaugeField_ptr()const{ return Deo_.getGaugeField_ptr(); }

  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;

  const Field md_force( const Field& eta,const Field& zeta) const;
  ////
  const Field mult_hop5_inv(const Field& f5) const;

  const Field gamma5(const Field& f5) const{ return Deo_.gamma5(f5);}
  ////
  const Field mult_eo(const Field& f) const; 
  const Field mult_oe(const Field& f) const; 

  const Field mult_eo_dag(const Field& f) const;
  const Field mult_oe_dag(const Field& f) const;

  const Field mult_oo(const Field& f)const;
  const Field mult_oo_inv(const Field& f)const;
  const Field mult_oo_dinv(const Field& f)const;

  const Field mult_ee(const Field& f)const;
  const Field mult_ee_inv(const Field& f)const;
  const Field mult_ee_dinv(const Field& f)const;

  void update_internal_state(){} 

};

#endif
