/*----------------------------------------------------------------------
!@filename fopr_signH_Zolotarev.h
Time-stamp: <2014-07-02 21:34:48 noaki>
------------------------------------------------------------------------*/
#ifndef FOPR_SIGN_ZOLOTAREV_INCLUDED
#define FOPR_SIGN_ZOLOTAREV_INCLUDED

#include <vector>
#include "Fopr/fopr.h"
#include "include/field.h"
#include "Solver/multiShiftSolver_CG.hpp"
#include "EigenModes/eigenModesSolver.hpp"
#include "EigenModes/eigenSorter.hpp"
#include "include/format_F.h"

struct EigenData{
  
  EigenData():Nlow(0){}
  // eigenmodes
  std::vector<double> lmd;
  std::vector<Field> evec;
  int Nlow;
  
  // Zolotarev coeffs
  std::vector<double> bl;
  std::vector<double> cl;
  std::vector<double> sigma;
};

struct EigenPrms{

  EigenPrms(int _Nmm,
	    int _Nk_l,
	    int _Nk_h,
	    int _Np_l,
	    int _Np_h,
	    int _Niter_l,
	    int _Niter_h,
	    double _enorm_l,
	    double _enorm_h,
	    double _vthrs_l,
	    double _vthrs_h,
	    double _stp_cnd,
	    int _Niter,
	    int _Npoly)
 :Nmm(_Nmm),Nk_l(_Nk_l),Nk_h(_Nk_h),Np_l(_Np_l),Np_h(_Np_h),
  Niter_l(_Niter_l),Niter_h(_Niter_h),
  enorm_l(_enorm_l),enorm_h(_enorm_h),
  vthrs_l(_vthrs_l),vthrs_h(_vthrs_h),
  stp_cnd(_stp_cnd),Niter(_Niter),Npoly(_Npoly){}
  
  // for calc of eigenmodes
  int Nmm;
  int Nk_l,Nk_h;
  int Np_l,Np_h;
  int Niter_l,Niter_h;
  double enorm_l,enorm_h;
  double vthrs_l,vthrs_h;

  // for Zolotarev approx
  double stp_cnd;
  int Niter;
  int Npoly;
};

class Fopr_signH_Zolotarev: public Fopr {
  const DiracWilsonLike* D_;
  const MultiShiftSolver_CG* msslv_;
  EigenData* const ed_;
  
  void evaluate_lowmodes(Field&, const Field&) const;
  void subtract_lowmodes(Field&) const;

  const Field multDdagD(const Field& f) const{return D_->mult_dag(D_->mult(f));}
  const int Np_;

public:
  Fopr_signH_Zolotarev(const DiracWilsonLike* D,
		       const EigenPrms& Eprms, 
		       EigenData* const ed)
    :ed_(ed),D_(D),Np_(Eprms.Npoly),
     msslv_(new MultiShiftSolver_CG(new Fopr_DdagD(D_),
				    Eprms.stp_cnd,Eprms.Niter)){}
  
  ~Fopr_signH_Zolotarev(){ delete msslv_;}
  const Field mult(const Field& f) const;
  const Field mult_dag(const Field& f) const{return mult(f);}
  const Field gamma5(const Field&f) const { return D_->gamma5(f); }

  const Field* getGaugeField_ptr()const{ return D_->getGaugeField_ptr();}

  void calc_force(Field& force,const Field& eta,const Field& zeta) const; 

  size_t fsize() const{return D_->fsize();}
  size_t gsize() const{return D_->gsize();}
};

#endif

