//---------------------------------------------------------
/*
  @file dirac_wilson_EvenOdd.h
  @brief Definition of Even Odd wilson operator
*/
//---------------------------------------------------------
#ifndef DIRAC_WILSON_EVENODD_INCLUDED
#define DIRAC_WILSON_EVENODD_INCLUDED

#include "dirac_wilson.h"
//#include "Communicator/comm_io.hpp"

namespace Dw{
  class EvOd:private Dirac_Wilson{
  private:
    SiteIndex_eo* idx_;
    int gauge_site_p(int hsite)const { return idx_->esec(hsite);}
    int gauge_site_m(int hsite)const { return idx_->osec(hsite);}
  public:    
    EvOd(double mass,const Field* u)
      :Dirac_Wilson(mass,u,EOtag()),
       idx_(SiteIndex_eo::instance()){}
    size_t fsize() const{return Dirac_Wilson::fsize();}
    size_t gsize() const{return Dirac_Wilson::gsize();}
    
    const Format::Format_F get_fermionFormat() const {
      return Dirac_Wilson::get_fermionFormat();
    }
    const Field mult(const Field&) const;
    const Field mult_dag(const Field&) const;

    //Preconditioning methods
    const Field mult_prec(const Field&)const{}//empty now preconditioned 
    const Field mult_dag_prec(const Field&)const{}//empty now preconditioned
    const Field left_precond(const Field&)const{}//empty now
    const Field right_precond(const Field&)const{}//empty now
    ////////////////////////////////////////////////////////////////////

    const Field gamma5(const Field& f) const{
      return Dirac_Wilson::gamma5(f);
    }
    const Field md_force(const Field& eta,const Field& zeta) const;
  };

  class OdEv:private Dirac_Wilson{
  private:
    SiteIndex_eo* idx_;
    int gauge_site_p(int hsite)const { return idx_->osec(hsite);}
    int gauge_site_m(int hsite)const { return idx_->esec(hsite);}
  public:    
    OdEv(double mass,const Field* u):Dirac_Wilson(mass,u,OEtag()),
				     idx_(SiteIndex_eo::instance()){}
    size_t fsize() const{return Dirac_Wilson::fsize();}
    size_t gsize() const{return Dirac_Wilson::gsize();}

    const Format::Format_F get_fermionFormat() const {
      return Dirac_Wilson::get_fermionFormat();
    }
    const Field mult(const Field&) const;
    const Field mult_dag(const Field&) const;

    //Preconditioning methods
    const Field mult_prec(const Field&)const{}//empty now preconditioned 
    const Field mult_dag_prec(const Field&)const{}//empty now preconditioned
    const Field left_precond(const Field&)const{}//empty now
    const Field right_precond(const Field&)const{}//empty now
    ///////////////////////////////////////////////////////////////////

    const Field gamma5(const Field& f) const{
      return Dirac_Wilson::gamma5(f);
    }
    const Field md_force(const Field& eta,const Field& zeta) const;
  };
}

class Dirac_Wilson_EvenOdd:public DiracWilsonLike_EvenOdd {
private:
  const Dw::EvOd Deo_;
  const Dw::OdEv Doe_;

public:
  Dirac_Wilson_EvenOdd(double mass,const Field* u)
    :Deo_(mass,u),Doe_(mass,u){}

  Dirac_Wilson_EvenOdd(const XML::node& node,const Field* u)
    :Deo_(Dw::read_mass(node),u),Doe_(Dw::read_mass(node),u){
    std::cout<<"Dirac_Wilson_EvenOdd created"<<std::endl;
  }
  
  size_t fsize() const{return Deo_.fsize();}
  size_t gsize() const{return Deo_.gsize();}

  const Field gamma5(const Field& f) const{return Deo_.gamma5(f);}

  const Field mult(const Field&) const;
  const Field mult_dag(const Field&) const;

  //Preconditioning methods
  const Field mult_prec(const Field&)const{}//empty now preconditioned 
  const Field mult_dag_prec(const Field&)const{}//empty now preconditioned
  const Field left_precond(const Field&)const{}//empty now
  const Field right_precond(const Field&)const{}//empty now
  //////////////////////////////////////////////////////////////////////
  

  const Field md_force(const Field&,const Field&) const;

  const Field mult_eo(const Field& f) const {return Deo_.mult(f);}
  const Field mult_oe(const Field& f) const {return Doe_.mult(f);}
  const Field mult_eo_dag(const Field& f) const {
    return Deo_.mult_dag(f);
  }
  const Field mult_oe_dag(const Field& f) const {
    return Doe_.mult_dag(f);
  }
  const Field mult_oo(const Field& f)const {return f;}
  const Field mult_ee(const Field& f)const {return f;}
  const Field mult_oo_inv(const Field& f)const {return f;}
  const Field mult_ee_inv(const Field& f)const {return f;}

  const ffmt_t get_fermionFormat() const {
    return Deo_.get_fermionFormat();
  }
  const Field operator()(int OpType,const Field&)const{}
};

#endif
