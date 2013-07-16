/*!
 * @file dirac_overlap_Zolotarev.h
 *
 * @brief Definition of Dirac_overlap_Zolotarev class
 */
#ifndef DIRAC_OVERLAP_Zolotarev_INCLUDED
#define DIRAC_OVERLAP_Zolotarev_INCLUDED

#include "include/fopr_signH_Zolotarev.h"

struct Dirac_overlap_ZolotarevParams{
  double M0;
  double mq;
  
  Dirac_overlap_ZolotarevParams(XML::node node) {
    XML::read(node, "M0", M0);
    XML::read(node, "mass", mq);
  }
  
  Dirac_overlap_ZolotarevParams(double M0_,double mq_)
    :M0(M0_),mq(mq_){}
};

class Dirac_overlap_Zolotarev : public DiracWilsonLike {
private:
  const Dirac_overlap_ZolotarevParams Params;
  const Fopr_signH_Zolotarev* Fopr_signH_;

public:
  Dirac_overlap_Zolotarev(XML::node node,
			  const Fopr_signH_Zolotarev* Fopr_signH)
    :Params(Dirac_overlap_ZolotarevParams(node)),
     Fopr_signH_(Fopr_signH){}

  Dirac_overlap_Zolotarev(double M0,double mass,
			  const Fopr_signH_Zolotarev* Fopr_signH)
    :Params(Dirac_overlap_ZolotarevParams(M0, mass)),
     Fopr_signH_(Fopr_signH){}
  
  size_t fsize()const {return Fopr_signH_->fsize();}
  size_t gsize()const {return Fopr_signH_->gsize();}
  double getMass()const{return Params.mq;}
  
  const Field* getGaugeField_ptr()const{ 
    return Fopr_signH_->getGaugeField_ptr(); }
  
  const Field mult    (const Field& f) const;
  const Field mult_dag(const Field& f) const;
  const Field gamma5  (const Field& f) const;

  const Field md_force(const Field& eta,const Field& zeta) const;
  void update_internal_state(){}
};

#endif
