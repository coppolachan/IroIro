/*! @file dirac_WilsonLike.hpp
 *  @brief defines the abstract base classes for DiracWilsonLike 
 * Time-stamp: <2013-05-23 10:08:18 noaki>
 */
#ifndef DIRAC_WILSONLIKE_INCLUDED
#define DIRAC_WILSONLIKE_INCLUDED

#include "dirac.hpp"
#include "include/format_F.h"
#include "include/format_G.h"

typedef Format::Format_F ffmt_t;
typedef Format::Format_G gfmt_t;

namespace Dop{
  // for EvenOdd
  struct EOtag{};    
  struct OEtag{};    

  double read_mass(const XML::node& node);
}

/*!
 * @class DiracWilsonLike
 * @brief Declaration of abstract base class for Dirac operators of Wilson type
 */
class DiracWilsonLike : public Dirac {

protected:
  //#ifdef IMPROVED_WILSON
  int r0(int c)const{return 2*c;}
  int r1(int c)const{return 2*(NC_+c);}
  int r2(int c)const{return 2*(2*NC_+c);}
  int r3(int c)const{return 2*(3*NC_+c);} 

  int i0(int c)const{return 2*c+1;}
  int i1(int c)const{return 2*(NC_+c)+1;}
  int i2(int c)const{return 2*(2*NC_+c)+1;}
  int i3(int c)const{return 2*(3*NC_+c)+1;} 

  int re(int c1,int c2)const{return 2*(NC_*c1+c2);}
  int im(int c1,int c2)const{return 2*(NC_*c1+c2)+1;}

  void gammaXcore(double*,const double*)const;
  void gammaYcore(double*,const double*)const;
  void gammaZcore(double*,const double*)const;
  void gammaTcore(double*,const double*)const;
  void gamma5core(double*,const double*)const;
  void projPcore(double*,const double*)const;
  void projMcore(double*,const double*)const;
  //#endif
public:
  virtual ~DiracWilsonLike(){}

  virtual const Field gamma5(const Field&) const = 0;
  virtual void md_force_p(Field&,const Field&,const Field&)const{}
  virtual void md_force_m(Field&,const Field&,const Field&)const{}

#ifdef IBM_BGQ_WILSON
  virtual void mult_ptr(double*,double* const)const{}
  virtual void mult_dag_ptr(double*,double* const)const{}  
  virtual void mult_ptr_EO(double*,double* const)const{}
  virtual void mult_dag_ptr_EO(double*,double* const)const{}  
#endif

}; 

/*!
 * @class DiracWilsonLike_EvenOdd
 * @brief Declaration of abstract base class for Dirac operators of Wilson type
 */
class DiracWilsonLike_EvenOdd : public DiracWilsonLike {
public:
  virtual ~DiracWilsonLike_EvenOdd(){}

  virtual const Field mult_eo(const Field&) const =0;
  virtual const Field mult_oe(const Field&) const =0;
  virtual const Field mult_eo_dag(const Field&) const =0;
  virtual const Field mult_oe_dag(const Field&) const =0;

  virtual const Field mult_ee(const Field&) const =0;
  virtual const Field mult_oo(const Field&) const =0;
  virtual const Field mult_ee_inv(const Field&) const =0;
  virtual const Field mult_oo_inv(const Field&) const =0;

  virtual const DiracWilsonLike* getDeo() const =0;
  virtual const DiracWilsonLike* getDoe() const =0;
};

/*!
 * @class Dirac_OptimalDomainWall_4D
 * @brief Declaration of abstract base class for 4D-reducted Domain-Wall fermions
 */
class Dirac_optimalDomainWall_4D : public DiracWilsonLike {
protected:
  void BprojCore(double* f,const double* f0,const double* fN)const;
  void BprojCore_dag(double* f0,double* fN,const double* f) const;
public:
  virtual ~Dirac_optimalDomainWall_4D(){}

  virtual double getMass()const =0;
  virtual const Field mult_inv(const Field&) const =0;
  virtual const Field mult_dag_inv(const Field&) const =0;
};

#endif
