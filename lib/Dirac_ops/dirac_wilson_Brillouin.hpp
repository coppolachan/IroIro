/*! @filename dirac_wilson_Brillouin.hpp
 * @brief declaration of Dirac_Wilson_Brillouin class 
 * Time-stamp: <2013-04-23 17:35:48 noaki>
 ----------------------------------------------------------------------*/
#ifndef DIRAC_WILSON_BRILLOUIN_INCLUDED
#define DIRAC_WILSON_BRILLOUIN_INCLUDED

#include "dirac_WilsonLike.hpp"
#include "include/pugi_interface.h"

#include "Main/Geometry/siteIndex.hpp"
#include "Main/Geometry/siteMap.hpp"

enum ImpType{Standard,Improved};

class Dirac_Wilson_Brillouin: public DiracWilsonLike{
private:
  int Nvol_,Nin_;
  double kbr_,m_;
  const Communicator* comm_;

  const Field* const u_;
  const ffmt_t ff_;
  const gfmt_t gf_;
  const size_t fsize_;

  int sr(int s,int c)const{return 2*(s*NC_+c);}
  int si(int s,int c)const{return 2*(s*NC_+c)+1;}

  // full site operation
  int xsl(int x,int n,int dir)const{return SiteMap::shiftSite.xslice(x,n,dir);}
  int slsize(int dir)const{return SiteMap::shiftSite.slice_size(0,dir);}

  const Field sft_p(const Field& f,int dir)const;
  const Field sft_m(const Field& f,int dir)const;
  const Field lap(const Field& f,int dir,double a)const;
  const Field del(const Field& f,int dir,double a)const;

  const Field gamma_x(const Field&)const;
  const Field gamma_y(const Field&)const;
  const Field gamma_z(const Field&)const;
  const Field gamma_t(const Field&)const;

  static const Field (Dirac_Wilson_Brillouin::*gm[])(const Field&)const;

  void(Dirac_Wilson_Brillouin::*mult_core)(Field&,const Field&)const;
  void mult_std(Field&,const Field&)const;
  void mult_imp(Field&,const Field&)const;
  
  int sgm(int,int,int)const;
public:
 /*! @brief constructor to create instance with normal site indexing */
  Dirac_Wilson_Brillouin(double mass,const Field* u,ImpType imp=Standard)
    :Nvol_(CommonPrms::instance()->Nvol()),
     mult_core(&Dirac_Wilson_Brillouin::mult_std),
     kbr_(1.0/(mass+15.0/8.0)),u_(u),
     comm_(Communicator::instance()),
     ff_(Nvol_),fsize_(ff_.size()),gf_(Nvol_),
     Nin_(ff_.Nin()){
    //
    if(imp==Improved){
      kbr_= 1.0/(mass+225.0/128.0);
      mult_core = &Dirac_Wilson_Brillouin::mult_imp;
    }
  }

  Dirac_Wilson_Brillouin(const XML::node& node,const Field* u,
			 ImpType imp=Standard)
    :Nvol_(CommonPrms::instance()->Nvol()),
     mult_core(NULL),u_(u),
     comm_(Communicator::instance()),
     ff_(Nvol_),fsize_(ff_.size()),gf_(Nvol_),
     Nin_(ff_.Nin()){
    //
    double mass;
    XML::read(node,"mass",mass,MANDATORY);

    if(      imp==Standard){
      mult_core = &Dirac_Wilson_Brillouin::mult_std;
      kbr_= 1.0/((15.0/8.0)+mass);
    }else if(imp==Improved){          
      mult_core = &Dirac_Wilson_Brillouin::mult_imp;
      kbr_= 1.0/(mass+225.0/128.0);
    }
    if(mult_core==NULL){
      std::cerr<<"Failed to costruct Dirac_Wilson_Brillouin\n";
      abort;
    }
  }
  size_t fsize() const{return fsize_;}
  size_t gsize() const{return gf_.size();}

  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;
  const Field mult_H(const Field&)const;
  const Field mult_del(const Field&)const;
  const Field mult_lap(const Field&)const;

  ////////////////////////////////////////Preconditioned versions
  // Wilson operator has no defined preconditioner now 
  const Field mult_prec     (const Field& f)const{return f;}
  const Field mult_dag_prec (const Field& f)const{return f;}
  const Field left_prec     (const Field& f)const{return f;}
  const Field right_prec    (const Field& f)const{return f;}
  const Field left_dag_prec (const Field& f)const{return f;}
  const Field right_dag_prec(const Field& f)const{return f;}
  //////////////////////////////////////////////////////////////

  const Field gamma5(const Field&) const;
  double getKappa() const {return kbr_;} 
  void update_internal_state(){}
};

#endif
