/*! @file dirac_wilson_adjoint.hpp
 * @brief Dirac_Wilson_Adjoint class 
 Time-stamp: <2013-11-28 00:22:34 noaki>
 */
#ifndef DIRAC_WILSON_ADJOINT_INCLUDED
#define DIRAC_WILSON_ADJOINT_INCLUDED

#include "dirac_wilson.hpp"

class Dirac_Wilson_Adjoint: public DiracWilsonLike{
private:
  int Nvol_;
  double kpp_;
  Dirac_Wilson Dw_;
  GammaMatrix dm_; 
  const Field* const u_;

  const ffmt_t ff_;
  const afmt_t af_;
  const gfmt_t gf_;
  
  int r(int c)const{return 2*c;}
  int i(int c)const{return 2*c+1;}

  double laf11r(double* fp)const{ return fp[r(2)] +1.0/sqrt(3.0)*fp[r(7)];}
  double laf11i(double* fp)const{ return fp[i(2)] +1.0/sqrt(3.0)*fp[i(7)];}
  
  double laf12r(double* fp)const{ return fp[r(0)] +fp[i(1)];}
  double laf12i(double* fp)const{ return fp[i(0)] -fp[r(1)];}
  
  double laf13r(double* fp)const{ return fp[r(3)] +fp[i(4)];}
  double laf13i(double* fp)const{ return fp[i(3)] -fp[r(4)];}
  
  double laf21r(double* fp)const{ return fp[r(0)] -fp[i(1)];}
  double laf21i(double* fp)const{ return fp[i(0)] +fp[r(1)];}

  double laf22r(double* fp)const{ return -fp[r(2)] +1.0/sqrt(3.0)*fp[r(7)];}
  double laf22i(double* fp)const{ return -fp[i(2)] +1.0/sqrt(3.0)*fp[i(7)];}
  
  double laf23r(double* fp)const{ return fp[r(5)] +fp[i(6)];}
  double laf23i(double* fp)const{ return fp[i(5)] -fp[r(6)];}

  double laf31r(double* fp)const{ return fp[r(3)] -fp[i(4)];}
  double laf31i(double* fp)const{ return fp[i(3)] +fp[r(4)];}

  double laf32r(double* fp)const{ return fp[r(5)] -fp[i(6)];}
  double laf32i(double* fp)const{ return fp[i(5)] +fp[r(6)];}

  double laf33r(double* fp)const{ return -2.0/sqrt(3.0)*fp[r(7)];}
  double laf33i(double* fp)const{ return -2.0/sqrt(3.0)*fp[i(7)];}

  static double (Dirac_Wilson_Adjoint::*laf[])(double* fp)const;
  const Field lmd_ad2fd(const Field& f,int c)const;
  const Field lmd_ad2Ufd(const Field& f,int mu,int c)const;

  void lfa(std::valarray<double>& h,double* fp,int c0)const;
  void lfUa(std::valarray<double>& h,double* fp,double* up,int c0)const;

  void lmd_fd2ad(Field& w,const Field& f,int c)const;
  void lmd_fd2Uad(Field& w,const Field& f,int mu,int c)const;

  void mult_p(Field& af,const Field& f,int dir)const;
  void mult_m(Field& af,const Field& f,int dir)const;

  void mult_full(Field&,const Field&)const;   /*! @brief  (1-kpp*D)*f */
  void mult_offdiag(Field&,const Field&)const;/*! @brief  -kpp*D*f */

  int gsite(int site)const {return site;}
  int esec(int hs)const {return SiteIndex_EvenOdd::instance()->esec(hs);}
  int osec(int hs)const {return SiteIndex_EvenOdd::instance()->osec(hs);}

  int(Dirac_Wilson_Adjoint::*gp)(int)const;
  int(Dirac_Wilson_Adjoint::*gm)(int)const;

  void(Dirac_Wilson_Adjoint::*mult_core)(Field&,const Field&)const;

  Dirac_Wilson_Adjoint(const Dirac_Wilson_Adjoint&);   /*!< @brief simple copy is prohibited.*/
  /*
  /// for test purpose
  const Field fadU(const Field&,int)const;
  const Field Ufad(const Field&,int)const;
  void test_unitarity(const Field&,int)const;
  */
public:
  /*! @brief constructor to create instance with e/o site indexing */
  Dirac_Wilson_Adjoint(double mass,const Field* u,Dop::EOtag)
    :Nvol_(CommonPrms::instance()->Nvol()/2),
     kpp_(0.5/(4.0+mass)),Dw_(mass,u,Dop::EOtag()),dm_(NADJ_),
     ff_(Nvol_),af_(Nvol_),gf_(Nvol_*2),u_(u),
     gp(&Dirac_Wilson_Adjoint::esec),gm(&Dirac_Wilson_Adjoint::osec),
     mult_core(&Dirac_Wilson_Adjoint::mult_offdiag){}

  Dirac_Wilson_Adjoint(double mass,const Field* u,Dop::OEtag)
    :Nvol_(CommonPrms::instance()->Nvol()/2),
     kpp_(0.5/(4.0+mass)),Dw_(mass,u,Dop::EOtag()),dm_(NADJ_),
     ff_(Nvol_),af_(Nvol_),gf_(Nvol_*2),u_(u),
     gp(&Dirac_Wilson_Adjoint::osec),gm(&Dirac_Wilson_Adjoint::esec),
     mult_core(&Dirac_Wilson_Adjoint::mult_offdiag){}

  /*! @brief constructor to create instance with normal site indexing */
  Dirac_Wilson_Adjoint(double mass,const Field* u)
    :Nvol_(CommonPrms::instance()->Nvol()),
     kpp_(0.5/(4.0+mass)),Dw_(mass,u),dm_(NADJ_),
     ff_(Nvol_),af_(Nvol_),gf_(Nvol_),u_(u),
     gp(&Dirac_Wilson_Adjoint::gsite),gm(&Dirac_Wilson_Adjoint::gsite),
     mult_core(&Dirac_Wilson_Adjoint::mult_full){}

  Dirac_Wilson_Adjoint(const XML::node& node,const Field* u)
    :Nvol_(CommonPrms::instance()->Nvol()),Dw_(0.0,u),dm_(NADJ_),
     ff_(Nvol_),af_(Nvol_),gf_(Nvol_),u_(u),
     gp(&Dirac_Wilson_Adjoint::gsite),gm(&Dirac_Wilson_Adjoint::gsite),
     mult_core(&Dirac_Wilson_Adjoint::mult_full){
    //
    double mass;
    XML::read(node, "mass", mass,MANDATORY);
    kpp_= 0.5/(4.0+mass);
  }
  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;

  const Field gamma5(const Field&) const;
  const Field md_force(const Field&,const Field&)const;
  void md_force_p(Field&,const Field&,const Field&)const;
  void md_force_m(Field&,const Field&,const Field&)const;

  double getKappa() const {return kpp_;}  
  double getMass() const {return 0.5/kpp_-4.0;}  

  size_t fsize()const{return af_.size();}
  size_t gsize()const{return gf_.size();}

  const Field* getGaugeField_ptr()const{ return u_; }

  void update_internal_state(){}
};


#endif
