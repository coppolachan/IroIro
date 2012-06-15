/*!
 * @file dirac_DomainWall_EvenOdd.hpp
 *
 * @brief Declaration of class Dirac_optimalDomainWall_EvenOdd (5d operator)
 */
#ifndef DIRAC_OPTIMALDOMAINWALL_EVENODD_INCLUDED
#define DIRAC_OPTIMALDOMAINWALL_EVENODD_INCLUDED

#include "dirac_DomainWall.hpp"

static  double mult_timer;
static  double multdag_timer;
/*!
 * @brief Defines the 5d Domain Wall operator with even/odd site indexing
 */
class Dirac_optimalDomainWall_EvenOdd : public DiracWilsonLike_EvenOdd {
  const Dirac_optimalDomainWall Deo_;
  const Dirac_optimalDomainWall Doe_;


  void md_force_eo(Field&,const Field&,const Field&)const;
  void md_force_oe(Field&,const Field&,const Field&)const;
  
  Dirac_optimalDomainWall_EvenOdd(const Dirac_optimalDomainWall_EvenOdd&);
  /*!< simple copy is prohibited */
public:
  Dirac_optimalDomainWall_EvenOdd(XML::node DWF_node,const Field* u,
				  DWFType Type=Standard)
    :Deo_(DWF_node,u,Dw::EOtag(),Type),
     Doe_(DWF_node,u,Dw::OEtag(),Type){
#if VERBOSITY>4
    CCIO::cout<<"Dirac_optimalDomainWall_Evenodd created"<<std::endl;
#endif
  }
  /*! @brief copy constractor to create Pauli-Villars operator */
  Dirac_optimalDomainWall_EvenOdd(const Dirac_optimalDomainWall_EvenOdd& D, 
				  DWFType Type=Standard)
    :Deo_(D.Deo_,Type), Doe_(D.Doe_,Type){}
  
  Dirac_optimalDomainWall_EvenOdd(double b,double c,double M0,double mq,
				  const std::vector<double>& omega,
				  const Field* u)
    :Deo_(b,c,M0,mq,omega,u,Dw::EOtag()),
     Doe_(b,c,M0,mq,omega,u,Dw::OEtag()){}
  
  ~Dirac_optimalDomainWall_EvenOdd(){
    CCIO::cout << "DWF Timer mult: "<< mult_timer << "\n";
    CCIO::cout << "DWF Timer multdag: "<< multdag_timer << "\n";
}
  
  size_t f4size()const{ return Deo_.f4size();}
  size_t fsize() const{ return Deo_.fsize(); }
  size_t gsize() const{ return Deo_.gsize(); }
  
  const Field gamma5(const Field& f5) const{ return Deo_.gamma5(f5);}
  const Field gamma5_4d(const Field& f4) const{ return Deo_.gamma5_4d(f4);}

  double getMass() const{return Deo_.getMass();}
  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;

  //Preconditioning methods
  const Field mult_prec    (const Field& f)const{return f;}//empty now
  const Field mult_dag_prec(const Field& f)const{return f;}//empty now
  const Field left_prec     (const Field& f)const{return f;}//empty now
  const Field left_dag_prec (const Field& f)const{return f;}//empty now
  const Field right_prec    (const Field& f)const{return f;}//empty now
  const Field right_dag_prec(const Field& f)const{return f;}//empty now
  //////////////////////////////////////////////////////////////////////

  const Field md_force( const Field& eta,const Field& zeta) const;


  ////
  const Field mult_hop5_inv(const Field& f5) const;

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

  const ffmt_t get_fermionFormat() const{return  Deo_.get_fermionFormat();}
  const std::vector<int> get_gsite() const { 
    return SiteIndex_EvenOdd::instance()->get_gsite();}

  void update_internal_state(){} 

  ////////////////////////////
#ifdef IBM_BGQ_WILSON
  void BGQ_EO_mult(Field&, const Field&) const;
#endif
};


#endif
