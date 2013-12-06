#ifndef DOMAINWALLUTILS_INCLUDED
#define DOMAINWALLUTILS_INCLUDED

#include <string>
#include <string.h>
#include "include/macros.hpp"
#include "include/pugi_interface.h"
#include "wilsonLikeUtils.hpp"
#include "include/format_F.h"
#include "include/field.h"

namespace DWF {

  struct EvenOdd_tag{
    int EOtag;
    EvenOdd_tag(int tag=0):EOtag(tag){}
  };

  const std::vector<double> getOmega(int Ns,double lmd_min,double lmd_max);
  double read_wilson_mass(const XML::node& node);
}

enum DWFType {Regular, PauliVillars};

/*!
 * @brief Container for parameter of the 5d Optimal Domain Wall operator
 * Parameters for \f[D_{\rm dwf}(m_q) =\omega D_W(-m_0)(1+cL(m_q))+(1-L(m_q))\f]
 */
struct DomainWallParams{
  std::vector<double> omega_;/*!< @brief Weights defining the approximation */
  double N5_;/*!< @brief the length in the 5th direction (must be even) */
  double b_; /*!< @brief scale factor (b!=1 for scaled Shamir H_T) */
  double c_; /*!< @brief the kernel (H_W (c=0) or H_T (c=1)) */
  double M0_;/*!< @brief wilson mass (must be negative) */
  double mq_;/*!< @brief Bare quark mass \f$m_q\f$ */
  std::vector<double> bs_, cs_;
  std::vector<double> dp_, dm_;
  std::vector<double> es_, fs_;

  void set_arrays();

  DomainWallParams(XML::node DWF_node,DWFType Type=Regular);
  DomainWallParams(const DomainWallParams& prms,DWFType Type=Regular);
  DomainWallParams(double b,double c,double M0,double mq,
		   const std::vector<double>& omega);
};

class Proj5D{
public:
  virtual void operator()(Field&,const Field&,int)const = 0;
  virtual ~Proj5D(){}
};

template <int NCOL=NC_> 
class ProjP: public Proj5D{
  int Nvol_,N5_;
  Format::Format_Fermion<NCOL> fmt_;
  GammaMatrix dm_;
public:  
  ProjP(int Nvol,int N5):Nvol_(Nvol),N5_(N5),fmt_(Nvol,N5),dm_(NCOL){}

  void operator()(Field& w,const Field& f5,int s)const{
    for(int site=0; site<Nvol_; ++site)
      dm_.projPcore(w.getaddr(fmt_.index(0,site)),
		    const_cast<Field&>(f5).getaddr(fmt_.index(0,site,s)));
  }
};

template <int NCOL=NC_> 
class ProjM: public Proj5D{
  int Nvol_,N5_;
  Format::Format_Fermion<NCOL> fmt_;
  GammaMatrix dm_;
public:  
  ProjM(int Nvol,int N5):Nvol_(Nvol),N5_(N5),fmt_(Nvol,N5),dm_(NCOL){}

  void operator()(Field& w,const Field& f5,int s)const{
    for(int site=0; site<Nvol_; ++site)
      dm_.projMcore(w.getaddr(fmt_.index(0,site)),
		    const_cast<Field&>(f5).getaddr(fmt_.index(0,site,s)));
  }
};

template <int NCOL=NC_> 
class Op5D{
  int Nvol_,N5_;
  Format::Format_Fermion<NCOL> fmt_;
  GammaMatrix dm_;

  const Field get4d(const Field& f5,int s) const{
    return Field(f5[fmt_.ex_slice(s)]); }
  
  void set5d(Field& f5,const Field& f4,int s) const{
    f5.set(fmt_.ex_slice(s),f4.getva()); }

public:
  Op5D(int Nvol,int N5):Nvol_(Nvol),N5_(N5),fmt_(Nvol,N5),dm_(NCOL){}
  
  const Field gamma5(const Field& f5) const{
    Field w5(f5.size()); 
    for(int s=0; s<N5_; ++s){
      for(int site=0; site<Nvol_; ++site)
	dm_.gamma5core(w5.getaddr(fmt_.index(0,site,s)),
		       const_cast<Field&>(f5).getaddr(fmt_.index(0,site,s)));
    }
    return w5; 
  }

  const Field R5(const Field& f5) const{
    Field w5(f5.size()); 
    for(int s=0; s<N5_; ++s) set5d(w5,get4d(f5,s),N5_-s-1);
    return w5; 
  }

  const Field R5g5(const Field& f5) const{
    Field w5(f5.size());
    for(int s=0; s<N5_; ++s){
      for(int site=0; site<Nvol_; ++site)
	dm_.gamma5core(w5.getaddr(fmt_.index(0,site,N5_-s-1)),
		       const_cast<Field&>(f5).getaddr(fmt_.index(0,site,s)));
    }
    return w5;
  }
};

typedef ProjP<NADJ_> ProjPadj;
typedef ProjM<NADJ_> ProjMadj;

typedef ProjP<NC_>   ProjPfdm;
typedef ProjM<NC_>   ProjMfdm;

typedef Op5D<NC_>   Op5Dfdm;
typedef Op5D<NADJ_> Op5Dadj;

#endif
