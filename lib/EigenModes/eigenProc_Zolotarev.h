/*!
 * @file eigenProc_Zolotarev.h
 *
 * @brief Declaration of EigenProc_Zolotarev class
 *
 */
#ifndef EIGENMONITOR_INCLUDED
#define EIGENMONITOR_INCLUDED

#ifndef EIGENMODES_IRL_INCLUDED
#include "eigenModes_IRL.h"
#endif

#ifndef FOPR_INCLUDED
#include "fopr.h"
#endif

#ifndef FORMAT_G_INCLUDED
#include "format_G.h"
#endif

#include <vector>
class Field;

class EigenProc_Zolotarev{
private:
  Fopr_H* H_;
  SortEigen_low* sort_low_;
  SortEigen_high* sort_high_;

  EigenModes_IRL* eigen_low_;
  EigenModes_IRL* eigen_high_;
  size_t fsize_;    
  int Nmm_;
  int Np_;

public:  
  EigenProc_Zolotarev(const DiracWilsonLike* D,
		      const EigenPrms& Eprms)
    :H_(new Fopr_H(D)),
     fsize_(H_->fsize()),
     sort_low_(new SortEigen_low),
     eigen_low_(new EigenModes_IRL(H_,sort_low_,
				   Eprms.Nk_l,Eprms.Np_l,Eprms.enorm_l,
				   Eprms.vthrs_l,Eprms.Niter_l)),
     sort_high_(new SortEigen_high),
     eigen_high_(new EigenModes_IRL(H_,sort_high_,
				    Eprms.Nk_h,Eprms.Np_h,Eprms.enorm_h,
				    Eprms.vthrs_h,Eprms.Niter_h)),
    Nmm_(Eprms.Nmm),Np_(Eprms.Npoly){}
  
  ~EigenProc_Zolotarev(){
    delete eigen_high_;
    delete eigen_low_;
    delete sort_high_;
    delete sort_low_;
    delete H_;
  }

  
  void calc(EigenData& ed);

  void Zolotarev_coeffs(std::vector<double>& bl,
			std::vector<double>& cl,
			std::vector<double>& sigma,
			int Np,double xmax,double xmin);
};

void EigenProc_Zolotarev::calc(EigenData& ed){

  int Nlow, Nhigh;

  ed.lmd.resize(Nmm_);
  ed.evec.resize(Nmm_);
  for(int k=0;k<Nmm_;++k) ed.evec[k].resize(fsize_);
  
  eigen_high_->calc(ed.lmd,ed.evec,Nhigh);
  double lmd_max = ed.lmd[0];
  
  eigen_low_->calc(ed.lmd,ed.evec,Nlow);
  double lmd_min = ed.lmd[0];
  
  Zolotarev_coeffs(ed.bl,ed.cl,ed.sigma,Np_,lmd_max,lmd_min);
}

double sign_Zolotarev(double x,const std::vector<double>& bl, 
		      const std::vector<double>& cl);

void poly_Zolotarev(std::vector<double>& bl,std::vector<double>& cl,
		    double bmax,double& UK);

void Jakobi_elliptic(double uu,double emmc,double &sn,double &cn,double &dn);

void EigenProc_Zolotarev::Zolotarev_coeffs(std::vector<double>& bl,
					   std::vector<double>& cl,
					   std::vector<double>& sigma,
					   int Np,double xmax,double xmin){
  double UK=0.0;

  bl.resize(Np_);
  cl.resize(2*Np_);

  poly_Zolotarev(bl,cl,xmax/xmin,UK);

  sigma.resize(Np);
  for(int i=0; i<Np; i++){
    printf(" %3d %12.4e %12.4e %12.4e\n",i, cl[i], cl[i+Np], bl[i]);
    sigma[i] = cl[2*i]*xmin*xmin;
  }
}

#endif
