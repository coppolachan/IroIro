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
  
  Dirac_overlap_ZolotarevParams(const XML::node node) {
    XML::read(node, "M0", M0);
    XML::read(node, "mass", mq);
  }
  
  Dirac_overlap_ZolotarevParams(const double M0_, 
				const double mq_):
    M0(M0_),mq(mq_){}
};

class Dirac_overlap_Zolotarev : public DiracWilsonLike {

private:
  const Dirac_overlap_ZolotarevParams Params;
  const Fopr_signH_Zolotarev* Fopr_signH;

public:
  Dirac_overlap_Zolotarev(const XML::node node,
			  const Fopr_signH_Zolotarev* Fopr_signH_)
    :Params(Dirac_overlap_ZolotarevParams(node)),
     Fopr_signH(Fopr_signH_){}

 Dirac_overlap_Zolotarev(const double M0, const double mass,
			  const Fopr_signH_Zolotarev* Fopr_signH_)
   :Params(Dirac_overlap_ZolotarevParams(M0, mass)),
    Fopr_signH(Fopr_signH_){}
  
  size_t fsize()const {return Fopr_signH->fsize();}
  size_t gsize()const {return Fopr_signH->gsize();}

  const Field operator()(int, const Field&) const{}//temporary

  const Field mult    (const Field& f) const;
  const Field mult_dag(const Field& f) const;
  const Field gamma5  (const Field& f) const;

  ////////////////////////////////////////Preconditioned versions
  // 4d operator has no preconditioner now 
  const Field mult_prec     (const Field&f)const{return f;}
  const Field mult_dag_prec (const Field&f)const{return f;}
  const Field left_prec     (const Field&f)const{return f;}
  const Field right_prec    (const Field&f)const{return f;}
  const Field left_dag_prec (const Field&f)const{return f;}
  const Field right_dag_prec(const Field&f)const{return f;}
  //////////////////////////////////////////////////////////////

  const Field md_force(const Field& eta,const Field& zeta) const;
  const Field md_force_core(const Field& eta,const Field& zeta) const{};

  const Format::Format_F get_fermionFormat() const{
    return Fopr_signH->get_fermionFormat();
  }
  const std::vector<int> get_gsite() const {
    return Fopr_signH->get_gsite();
  }

  void update_internal_state(){};
};

#endif
