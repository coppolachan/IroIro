#include "Fopr/fopr_signH_Zolotarev.h"
#include "Fields/field_expressions.hpp"
#include <cassert>
#include <cmath>

const Field Fopr_signH_Zolotarev::mult(const Field& f) const{
  using namespace FieldExpression;

  int Nsbt = ed_->Nlow;
  printf("  sign function: Nsbt = %d\n",Nsbt);

  // Low-mode subtraction
  Field f2(f);
  if(Nsbt> 0) subtract_lowmodes(f2);

  //multi-shift solver
  std::vector<Field> xq(Np_);
  for(int i=0; i<Np_; ++i) xq[i].resize(D_->fsize());

  int Nconv;
  double diff;
  printf ("Shift solver\n");
  printf ("  Number of shift values = %d\n",Np_);

  msslv_->solve(xq,f2,ed_->sigma,diff,Nconv);
  
  Field v(D_->fsize());
  for(int i=0; i<Np_; ++i) v += ed_->bl[i]*xq[i];

  double lmd_min = ed_->lmd[0];

  Field w = multDdagD(v);
  w += ed_->cl[2*Np_-1]*lmd_min*lmd_min*v;

  v = 1.0/lmd_min*D_->gamma5(D_->mult(w));
  if(Nsbt> 0) evaluate_lowmodes(v,f);
  return v;
}

void Fopr_signH_Zolotarev::subtract_lowmodes(Field& w) const{
  assert(w.size()%2==0);

  double prdr, prdi;
  int Nsbt = ed_->Nlow;

  for(int k=0; k< Nsbt; ++k){
    prdr = 0.0;
    prdi = 0.0;

    for(int i=0; i< w.size(); i+=2){
      prdr += ed_->evec[k][i]*w[i  ] +ed_->evec[k][i+1]*w[i+1];
      prdi += ed_->evec[k][i]*w[i+1] -ed_->evec[k][i+1]*w[i  ];
    }
    for(int i=0; i< w.size(); i+=2){
      w.add(i  ,-prdr*ed_->evec[k][i  ] +prdi*ed_->evec[k][i+1]);
      w.add(i+1,-prdr*ed_->evec[k][i+1] -prdi*ed_->evec[k][i  ]);
    }
  }
}

void Fopr_signH_Zolotarev::
evaluate_lowmodes(Field& x,const Field& w) const{
  assert(w.size()%2==0);

  double prdr, prdi;
  int Nsbt = ed_->Nlow;

  for(int k=0; k< Nsbt; ++k){
    prdr = 0.0;
    prdi = 0.0;

    for(int i=0; i< w.size(); i+=2){
      prdr += ed_->evec[k][i]*w[i  ] +ed_->evec[k][i+1]*w[i+1];
      prdi += ed_->evec[k][i]*w[i+1] -ed_->evec[k][i+1]*w[i  ];
    }
    double sgn = ed_->lmd[k]/fabs(ed_->lmd[k]);
    prdr *= sgn;
    prdi *= sgn;
    
    for(int i=0; i< w.size(); i+=2){
      x.add(i  ,-prdr*ed_->evec[k][i  ] +prdi*ed_->evec[k][i+1]);
      x.add(i+1,-prdr*ed_->evec[k][i+1] -prdi*ed_->evec[k][i  ]);
    }
  }
}

void Fopr_signH_Zolotarev::
calc_force(Field& fce,const Field& eta,const Field& zeta) const{ 

  using namespace FieldExpression;
  int Nconv;
  double diff;
  double lmd_min = ed_->lmd[0];

  std::vector<Field> etq(Np_);
  for(int i=0; i<Np_; ++i) etq[i].resize(D_->fsize());
  msslv_->solve(etq,eta,ed_->sigma,diff,Nconv);
  for(int i=0; i<Np_; ++i) etq[i]*= lmd_min*lmd_min;

  std::vector<Field> ztq(Np_);
  for(int i=0; i<Np_; ++i) ztq[i].resize(D_->fsize());
  msslv_->solve(ztq,zeta,ed_->sigma,diff,Nconv);
  for(int i=0; i<Np_; ++i) ztq[i]*= lmd_min*lmd_min;

  double Bl =0.0;
  for(int i=0; i<Np_; ++i) Bl+=ed_->bl[i];

  fce = D_->md_force(eta,zeta);
  fce*= Bl/lmd_min;

  for(int i=0; i<Np_; ++i){
    double cl2i = ed_->cl[2*i];
    double coeff = ed_->bl[i]*cl2i*(ed_->cl[2*Np_-1]-cl2i)/lmd_min;
    fce+= coeff*D_->md_force(etq[i],ztq[i]);
  }
  for(int i=0; i<Np_; ++i){
    double coeff = ed_->bl[i]*(ed_->cl[2*Np_-1]-ed_->cl[2*i])/lmd_min;
    fce-= coeff*D_->md_force(D_->gamma5(D_->mult(etq[i])),
			     D_->gamma5(D_->mult(ztq[i])));
  }
  fce/= lmd_min*lmd_min;
}
