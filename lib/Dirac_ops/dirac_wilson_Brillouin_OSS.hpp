/*! @filename dirac_wilson_Brillouin_OSS.hpp
 * @brief declaration of Dirac_Wilson_Brillouin class 
 * Time-stamp: <2014-01-15 18:36:10 noaki>
 ----------------------------------------------------------------------*/
#ifndef DIRAC_WILSON_BRILLOUIN_OSS_INCLUDED
#define DIRAC_WILSON_BRILLOUIN_OSS_INCLUDED

#include "dirac_WilsonLike.hpp"
#include "include/pugi_interface.h"
#include "include/common_fields.hpp"
#include "Main/gaugeConf.hpp"
#include "wilsonLikeUtils.hpp"
//#include "Geometry/siteIndex.hpp"
#include "Geometry/siteMap.hpp"

class Dirac_Wilson_Brillouin_OSS: public DiracWilsonLike{
private:
  int Nvol_,Nin_;
  double kbr_;
  double mass_;
  const Communicator* comm_;

  const Field* const u_;
  const ffmt_t ff_;
  const gfmt_t gf_;
  const size_t fsize_;
  const gfmt_t Wfmt_;

  GammaMatrix dm_;
  Field W_;
  mutable std::vector<Field> F_;

  ///
  int sr(int s,int c)const{return 2*(s*NC_+c);}
  int si(int s,int c)const{return 2*(s*NC_+c)+1;}

  const Field der_iso()const;
  const Field der_iso_X()const;
  const Field der_iso_Y()const;
  const Field der_iso_Z()const;
  const Field der_iso_T()const;
  const Field lap_bri()const;
  void mult_col(FermionField&,int)const;

  static const Field (Dirac_Wilson_Brillouin_OSS::*der[])()const;
  static void (Dirac_Wilson_Brillouin_OSS::*ghops[])(GaugeField1D&,
						     const GaugeField&)const;

  void(Dirac_Wilson_Brillouin_OSS::*mult_core)(Field&,const Field&)const;
  void mult_std(Field&,const Field&)const;
  void mult_imp(Field&,const Field&)const;
  
  int sgm(int,int,int)const;
  void Wsetup(const GaugeField&);
  void Fsetup(const FermionField&)const;

  void shift_fwd(FermionField&,const FermionField&,int)const;  
  void shift_bwd(FermionField&,const FermionField&,int)const;  

#ifdef IBM_BGQ_WILSON
  Field unit_g_;
  double* gptr_;
#endif
public:
 /*! @brief constructor to create instance with normal site indexing */
  Dirac_Wilson_Brillouin_OSS(double mass,const Field* u,ImpType imp=Standard)
    :Nvol_(CommonPrms::instance()->Nvol()),
     mult_core(&Dirac_Wilson_Brillouin_OSS::mult_std),
     mass_(mass),
     kbr_(1.0/(mass+15.0/8.0)),u_(u),
     comm_(Communicator::instance()),
     ff_(Nvol_),fsize_(ff_.size()),gf_(Nvol_),
     Nin_(ff_.Nin()),Wfmt_(Nvol_,80),
#ifdef IBM_BGQ_WILSON     
     unit_g_(gf_.size()),gptr_(unit_g_.getaddr(0)),
#endif
     W_(Wfmt_.size()),F_(Wfmt_.Nex()){
    GaugeField gu(*u,Nvol_);
    Wsetup(gu);
    for(int i = 0; i<F_.size(); ++i) F_[i].resize(fsize_);
#ifdef IBM_BGQ_WILSON     
    GaugeConf_unit gunit(gf_);
    gunit.init_conf(unit_g_);
#endif
    //
    if(imp==Improved){
      kbr_= 1.0/(mass+225.0/128.0);
      mult_core = &Dirac_Wilson_Brillouin_OSS::mult_imp;
    }
  }

  Dirac_Wilson_Brillouin_OSS(const XML::node& node,const Field* u,
			 ImpType imp=Standard)
    :Nvol_(CommonPrms::instance()->Nvol()),
     mult_core(NULL),u_(u),
     comm_(Communicator::instance()),
     ff_(Nvol_),fsize_(ff_.size()),gf_(Nvol_),
     Nin_(ff_.Nin()),Wfmt_(Nvol_,80),
#ifdef IBM_BGQ_WILSON     
     unit_g_(gf_.size()),gptr_(unit_g_.getaddr(0)),
#endif
     W_(Wfmt_.size()),F_(Wfmt_.Nex()){
    //
    GaugeField gu(*u,Nvol_);
    Wsetup(gu);
    for(int i = 0; i<F_.size(); ++i) F_[i].resize(fsize_);
#ifdef IBM_BGQ_WILSON     
    GaugeConf_unit gunit(gf_);
    gunit.init_conf(unit_g_);
#endif
    //
    double mass;
    XML::read(node,"mass",mass,MANDATORY);

    if(      imp==Standard){
      mult_core = &Dirac_Wilson_Brillouin_OSS::mult_std;
      kbr_= 1.0/(mass+15.0/8.0);
    }else if(imp==Improved){          
      mult_core = &Dirac_Wilson_Brillouin_OSS::mult_imp;
      kbr_= 1.0/(mass+225.0/128.0);
    }
    if(mult_core==NULL){
      std::cerr<<"Failed to costruct Dirac_Wilson_Brillouin_OSS\n";
      abort;
    }
  }

  size_t fsize() const{return fsize_;}
  size_t gsize() const{return gf_.size();}
  const Field* getGaugeField_ptr()const{return u_; }
  
  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;

  const Field gamma(const Field&,int) const;
  const Field gamma5(const Field&) const;

  double getKappa() const {return kbr_;} 
  double getMass() const {return mass_;} 

private:

 //hopping terms
  //1-hop * 8
  void gpX(GaugeField1D&,const GaugeField&) const;
  void gmX(GaugeField1D&,const GaugeField&) const;
  void gpY(GaugeField1D&,const GaugeField&) const;
  void gmY(GaugeField1D&,const GaugeField&) const;
  void gpZ(GaugeField1D&,const GaugeField&) const;
  void gmZ(GaugeField1D&,const GaugeField&) const;
  void gpT(GaugeField1D&,const GaugeField&) const;
  void gmT(GaugeField1D&,const GaugeField&) const;
  //2-hop * 24
  void gpXpY(GaugeField1D&,const GaugeField&) const;
  void gpXpZ(GaugeField1D&,const GaugeField&) const;
  void gpXpT(GaugeField1D&,const GaugeField&) const;
  void gpYpZ(GaugeField1D&,const GaugeField&) const;
  void gpYpT(GaugeField1D&,const GaugeField&) const;
  void gpZpT(GaugeField1D&,const GaugeField&) const;
  void gmXmY(GaugeField1D&,const GaugeField&) const;
  void gmXmZ(GaugeField1D&,const GaugeField&) const;
  void gmXmT(GaugeField1D&,const GaugeField&) const;
  void gmYmZ(GaugeField1D&,const GaugeField&) const;
  void gmYmT(GaugeField1D&,const GaugeField&) const;
  void gmZmT(GaugeField1D&,const GaugeField&) const;
  void gpXmY(GaugeField1D&,const GaugeField&) const;
  void gpXmZ(GaugeField1D&,const GaugeField&) const;
  void gpXmT(GaugeField1D&,const GaugeField&) const;
  void gpYmZ(GaugeField1D&,const GaugeField&) const;
  void gpYmT(GaugeField1D&,const GaugeField&) const;
  void gpZmT(GaugeField1D&,const GaugeField&) const;
  void gmXpY(GaugeField1D&,const GaugeField&) const;
  void gmXpZ(GaugeField1D&,const GaugeField&) const;
  void gmXpT(GaugeField1D&,const GaugeField&) const;
  void gmYpZ(GaugeField1D&,const GaugeField&) const;
  void gmYpT(GaugeField1D&,const GaugeField&) const;
  void gmZpT(GaugeField1D&,const GaugeField&) const;

  //4-hop * 32
  void gpXpYpZ(GaugeField1D&,const GaugeField&) const;
  void gpXpYpT(GaugeField1D&,const GaugeField&) const;
  void gpXpZpT(GaugeField1D&,const GaugeField&) const;
  void gpYpZpT(GaugeField1D&,const GaugeField&) const;
  void gpXpYmZ(GaugeField1D&,const GaugeField&) const;
  void gpXpYmT(GaugeField1D&,const GaugeField&) const;
  void gpXpZmT(GaugeField1D&,const GaugeField&) const;
  void gpYpZmT(GaugeField1D&,const GaugeField&) const;
  void gpXmYpZ(GaugeField1D&,const GaugeField&) const;
  void gpXmYpT(GaugeField1D&,const GaugeField&) const;
  void gpXmZpT(GaugeField1D&,const GaugeField&) const;
  void gpYmZpT(GaugeField1D&,const GaugeField&) const;
  void gmXpYpZ(GaugeField1D&,const GaugeField&) const;
  void gmXpYpT(GaugeField1D&,const GaugeField&) const;
  void gmXpZpT(GaugeField1D&,const GaugeField&) const;
  void gmYpZpT(GaugeField1D&,const GaugeField&) const;
  void gpXmYmZ(GaugeField1D&,const GaugeField&) const;
  void gpXmYmT(GaugeField1D&,const GaugeField&) const;
  void gpXmZmT(GaugeField1D&,const GaugeField&) const;
  void gpYmZmT(GaugeField1D&,const GaugeField&) const;
  void gmXpYmZ(GaugeField1D&,const GaugeField&) const;
  void gmXpYmT(GaugeField1D&,const GaugeField&) const;
  void gmXpZmT(GaugeField1D&,const GaugeField&) const;
  void gmYpZmT(GaugeField1D&,const GaugeField&) const;
  void gmXmYpZ(GaugeField1D&,const GaugeField&) const;
  void gmXmYpT(GaugeField1D&,const GaugeField&) const;
  void gmXmZpT(GaugeField1D&,const GaugeField&) const;
  void gmYmZpT(GaugeField1D&,const GaugeField&) const;
  void gmXmYmZ(GaugeField1D&,const GaugeField&) const;
  void gmXmYmT(GaugeField1D&,const GaugeField&) const;
  void gmXmZmT(GaugeField1D&,const GaugeField&) const;
  void gmYmZmT(GaugeField1D&,const GaugeField&) const;

  //4-hop * 12
  void gpXpYpZpT(GaugeField1D&,const GaugeField&) const;
  void gpXpYpZmT(GaugeField1D&,const GaugeField&) const;
  void gpXpYmZpT(GaugeField1D&,const GaugeField&) const;
  void gpXmYpZpT(GaugeField1D&,const GaugeField&) const;
  void gmXpYpZpT(GaugeField1D&,const GaugeField&) const;
  void gpXpYmZmT(GaugeField1D&,const GaugeField&) const;
  void gpXmYpZmT(GaugeField1D&,const GaugeField&) const;
  void gpXmYmZpT(GaugeField1D&,const GaugeField&) const;
  void gmXpYpZmT(GaugeField1D&,const GaugeField&) const;
  void gmXpYmZpT(GaugeField1D&,const GaugeField&) const;
  void gmXmYpZpT(GaugeField1D&,const GaugeField&) const;
  void gpXmYmZmT(GaugeField1D&,const GaugeField&) const;
  void gmXpYmZmT(GaugeField1D&,const GaugeField&) const;
  void gmXmYpZmT(GaugeField1D&,const GaugeField&) const;
  void gmXmYmZpT(GaugeField1D&,const GaugeField&) const;
  void gmXmYmZmT(GaugeField1D&,const GaugeField&) const;
};

#endif
