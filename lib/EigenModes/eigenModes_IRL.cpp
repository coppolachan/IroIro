//---------------------------------------------------------------------
// eigenModes_IRL.cpp
//---------------------------------------------------------------------
#include "eigenModes_IRL.h"
#include "sortEigen.h"
#include "include/fopr.h"
#include <cassert>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "Communicator/comm_io.hpp"

using namespace std;

class Field;

void EigenModes_IRL::calc(vector<double>& lmd, vector<Field>& evec,
			  const Field& b, int& Nsbt, int& Nconv)const{
  using namespace FieldExpression;
  size_t fsize = evec[0].size();
  
  Nconv = -1;
  Nsbt = 0;
  int Nm = Nk_+Np_;

  CCIO::cout << "Nk = " << Nk_ << " Np = "<< Np_ << endl;
  CCIO::cout << "Nm = " << Nm << endl;
  CCIO::cout << "size of lmd   = " << lmd.size() << endl;
  CCIO::cout << "size of evec  = " << evec.size() << endl;
  
  assert(Nm < evec.size() && Nm < lmd.size());
  
  vector<double> lme(Nm);
  vector<double> lmd2(Nm);
  vector<double> lme2(Nm);
  vector<double> Qt(Nm*Nm);
  vector<int>    Iconv(Nm);
  
  vector<Field>  B(Nm);
  for(int k=0; k<Nm; ++k) B[k].resize(fsize);
  
  Field f(fsize);
  Field v(fsize);
  

  int k1 = 1;
  int k2 = Nk_;
  int kconv = 0;
  
  int Kdis  = 0;
  int Kthrs = 0;
  double beta_k;
  
  // Set initial vector
  evec[0] = 1.0;
  double vnorm = evec[0]*evec[0];
  evec[0] = 1.0/sqrt(vnorm);
  // (uniform vector)


  // Initial Nk steps
  for(int k=0; k<k2; ++k) step(lmd,lme,evec,f,Nm,k);

  // Restarting loop begins
  for(int iter = 0; iter<Niter_; ++iter){
    CCIO::cout<<"\n iter= "<< iter << endl;

    int Nm2 = Nm - kconv;
    for(int k=k2; k<Nm; ++k) step(lmd,lme,evec,f,Nm,k);
    f *= lme[Nm-1];

    // getting eigenvalues
    for(int k=0; k<Nm2; ++k){
      lmd2[k] = lmd[k+k1-1];
      lme2[k] = lme[k+k1-1];
    }
    setUnit_Qt(Nm,Qt);
    diagonalize(lmd2,lme2,Nm2,Nm,Qt);
    // sorting
    sort_->push(lmd2,Nm);

    // Implicitly shifted QR transformations
    setUnit_Qt(Nm,Qt);
    for(int ip=k2; ip<Nm; ++ip) 
      qr_decomp(lmd,lme,Nm,Nm,Qt,lmd2[ip],1,Nm);
    
    for(int i=0; i<(Nk_+1); ++i) B[i] = 0.0;

    for(int j=k1-1; j<k2+1; ++j){
      for(int k=0; k<Nm; ++k){
	B[j] += Qt[k+Nm*j] * evec[k];
      }
    }
    for(int j=k1-1; j<k2+1; ++j) evec[j] = B[j];

    // Compressed vector f and beta(k2)
    f *= Qt[Nm-1+Nm*(k2-1)];
    f += lme[k2-1] * evec[k2];
    beta_k = f * f;
    beta_k = sqrt(beta_k);
    CCIO::cout<<" beta(k) = "<<beta_k<<endl;

    double betar = 1.0/beta_k;
    evec[k2] = betar * f;
    lme[k2-1] = beta_k;

    // Convergence test
    for(int k=0; k<Nm2; ++k){    
      lmd2[k] = lmd[k];
      lme2[k] = lme[k];
    }
    setUnit_Qt(Nm,Qt);
    diagonalize(lmd2,lme2,Nk_,Nm,Qt);

    for(int k = 0; k<Nk_; ++k) B[k]=0.0;

    for(int j = 0; j<Nk_; ++j){
      for(int k = 0; k<Nk_; ++k){
	B[j] += Qt[k+j*Nm] * evec[k];
      }
    }
    Kdis = 0;
    Kthrs = 0;

    for(int i=0; i<Nk_; ++i){
      v = opr_->mult(B[i]);
      double vnum = B[i] * v;
      double vden = B[i] * B[i];
      lmd2[i] = vnum/vden;
      v -= lmd2[i] * B[i];
      double vv = v * v;
      CCIO::cout<<"  ["<<i<<"]  "<<lmd2[i]<<"  "<<vv<< endl;
      if(vv<enorm_){
        Iconv[Kdis] = i;
        ++Kdis;
        if(sort_->saturated(lmd2[i],vthrs_)) ++Kthrs;
      }
    }  // i-loop end

    CCIO::cout<<" #modes converged: "<<Kdis<<endl;

    if(Kthrs > 0){
      // (there is a converged eigenvalue larger than Vthrs.)
      Nconv = iter;
      goto converged;
    }
  } // end of iter loop

  CCIO::cout<<"\n NOT converged.\n";
  abort();

converged:
  // Sorting

  lmd.clear();
  evec.clear();
  for(int i=0; i<Kdis; ++i){
    lmd.push_back(lmd2[Iconv[i]]);
    evec.push_back(B[Iconv[i]]);
  }
  sort_->push(lmd,evec,Kdis);
  Nsbt = Kdis - Kthrs;

  printf("\n Converged:\n");
  printf("  Nconv   = %d\n",Nconv);
  printf("  beta(k) = %20.14e\n",beta_k);
  printf("  Kdis    = %d\n",Kdis);
  printf("  Nsbt    = %d\n",Nsbt);
}

/* ====================================================== */
void EigenModes_IRL::
step(vector<double>& lmd,
     vector<double>& lme,
     vector<Field>& evec, Field& w,int Nm, int k)const{

  using namespace FieldExpression;

  if(k>=Nm){
    // printf("k larger than Nm: irlegal.\n");
    abort();
  }else if(k==0){  // Initial step
    w = opr_->mult(evec[k]);

    double alph = evec[k] * w;
    w -= alph * evec[k];
    lmd[k] = alph;

    double beta = w * w;
    beta = sqrt(beta);
    double betar = 1.0/beta;
    evec[k+1] = betar * w;
    lme[k] = beta;

  }else{   // Iteration step
    w = opr_->mult(evec[k]);
    w -= lme[k-1] * evec[k-1];

    double alph = evec[k] * w;
    w -= alph * evec[k];

    double beta = w * w;
    beta = sqrt(beta);
    double betar = 1.0/beta;
    w *= betar;

    lmd[k] = alph;
    lme[k] = beta;
    orthogonalize(w,evec,k);

    if(k < Nm-1) evec[k+1] = w;
  }
}

/*
void EigenModes_IRL::orthogonalize(Field& w,const vector<Field>& evec,int k){
  // Schmidt orthogonalization                                   
                                                
  size_t size = w.size();
  assert(size%2 ==0);

  std::slice re(0,size/2,2);
  std::slice im(1,size/2,2);

  for(int j=0; j<k; ++j){
    double prdr = evec[j]*w;
    double prdi = evec[j].im_prod(w);

    valarray<double> evr(evec[j][re]);
    valarray<double> evi(evec[j][im]);

    w.add(re, -prdr*evr +prdi*evi);
    w.add(im, -prdr*evi -prdi*evr);
  }
}
*/

void EigenModes_IRL::
orthogonalize(Field& w,const vector<Field>& vk, int k)const{
  // Schmidt orthogonalization

  if( (w.size()%2)!=0 ) abort();

  double prdr, prdi;
  double vr, vi;

  for(int j=0; j<k; ++j){
    prdr = 0.0;
    prdi = 0.0;
    for(int i=0; i<w.size(); i+=2){
      prdr += vk[j][i]*w[i]  +vk[j][i+1]*w[i+1];
      prdi += vk[j][i]*w[i+1]-vk[j][i+1]*w[i];
    }

    for(int i=0; i<w.size(); i+=2){

      vr = w[i]  -prdr*vk[j][i]  +prdi*vk[j][i+1];
      vi = w[i+1]-prdr*vk[j][i+1]-prdi*vk[j][i];
      w.set(i  ,vr);
      w.set(i+1,vi);

      //vr = -prdr*vk[j][i]  +prdi*vk[j][i+1];
      //vi = -prdr*vk[j][i+1]-prdi*vk[j][i];
      //w.add(i  ,vr);
      //w.add(i+1,vi);
    }
  }
}



void EigenModes_IRL::
setUnit_Qt(int Nm, vector<double>& Qt)const{
  for(int i=0; i<Qt.size(); ++i) Qt[i] = 0.0;
  for(int k=0; k<Nm; ++k) Qt[k + k*Nm] = 1.0;
}

void EigenModes_IRL::
diagonalize(vector<double>& lmd, vector<double>& lme,
            int Nk, int Nm, vector<double>& Qt)const{

  int Niter = 100*Nm;
  int kmin = 1;
  int kmax = Nk;
  // (this should be more sophisticated)

  for(int iter=0; iter<Niter; ++iter){

    // determination of 2x2 leading submatrix
    double dsub = lmd[kmax-1]-lmd[kmax-2];
    double dd = sqrt(dsub*dsub + 4.0*lme[kmax-2]*lme[kmax-2]);
    double Dsh = 0.5*(lmd[kmax-2]+lmd[kmax-1] +dd*(dsub/fabs(dsub)));
    // (Dsh: shift)

    // transformation
    qr_decomp(lmd,lme,Nk,Nm,Qt,Dsh,kmin,kmax);

    // Convergence criterion (redef of kmin and kamx)
    for(int j=kmax-1; j>= kmin; --j){
      double dds = fabs(lmd[j-1])+fabs(lmd[j]);
      if(fabs(lme[j-1])+dds > dds){
        kmax = j+1;
        goto continued;
      }
    }
    Niter = iter;
    //    printf("  converged at iter = %d\n",Niter);
    return;

  continued:
    for(int j=0; j<kmax-1; ++j){
      double dds = fabs(lmd[j])+fabs(lmd[j+1]);
      if(fabs(lme[j])+dds > dds){
        kmin = j+1;
        break;
      }
    }
  }
  printf("QL method too many iteration: %d\n",Niter);
  abort();
}

void EigenModes_IRL::
qr_decomp(vector<double>& lmd, vector<double>& lme,
	  int Nk, int Nm, vector<double>& Qt,
	  double Dsh, int kmin, int kmax)const{
  int k = kmin-1;
  double x;

  double Fden = 1.0/sqrt((lmd[k]-Dsh)*(lmd[k]-Dsh) +lme[k]*lme[k]);
  double c = ( lmd[k] -Dsh) *Fden;
  double s = -lme[k] *Fden;

  double tmpa1 = lmd[k];
  double tmpa2 = lmd[k+1];
  double tmpb  = lme[k];

  lmd[k]   = c*c*tmpa1 +s*s*tmpa2 -2.0*c*s*tmpb;
  lmd[k+1] = s*s*tmpa1 +c*c*tmpa2 +2.0*c*s*tmpb;
  lme[k]   = c*s*(tmpa1-tmpa2) +(c*c-s*s)*tmpb;
  x        = -s*lme[k+1];
  lme[k+1] = c*lme[k+1];

  for(int i=0; i<Nk; ++i){
    double Qtmp1 = Qt[i+Nm*k  ];
    double Qtmp2 = Qt[i+Nm*(k+1)];
    Qt[i+Nm*k    ] = c*Qtmp1 - s*Qtmp2;
    Qt[i+Nm*(k+1)] = s*Qtmp1 + c*Qtmp2; 
  }

  // Givens transformations
  for(int k = kmin; k < kmax-1; ++k){
    double Fden = 1.0/sqrt( x*x +lme[k-1]*lme[k-1]);
    double c = lme[k-1]*Fden;
    double s = - x*Fden;

    double tmpa1 = lmd[k];
    double tmpa2 = lmd[k+1];
    double tmpb  = lme[k];


    lmd[k]   = c*c*tmpa1 +s*s*tmpa2 -2.0*c*s*tmpb;
    lmd[k+1] = s*s*tmpa1 +c*c*tmpa2 +2.0*c*s*tmpb;
    lme[k]   = c*s*(tmpa1-tmpa2) +(c*c-s*s)*tmpb;
    lme[k-1] = c*lme[k-1] -s*x;

    if(k != kmax-2){
      x = -s*lme[k+1];
      lme[k+1] = c*lme[k+1];
    }

    for(int i=0; i<Nk; ++i){
      double Qtmp1 = Qt[i+Nm*k    ];
      double Qtmp2 = Qt[i+Nm*(k+1)];
      Qt[i+Nm*k    ] = c*Qtmp1 -s*Qtmp2;
      Qt[i+Nm*(k+1)] = s*Qtmp1 +c*Qtmp2;
    }
  }
}

