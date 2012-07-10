//---------------------------------------------------------------------
// eigenModes_IRL.cpp
//---------------------------------------------------------------------
#include "eigenModes_IRL.hpp"
#include "sortEigen.h"
#include "include/fopr.h"
#include "Fields/field_expressions.hpp"
#include <cassert>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

void EigenModes_IRL::calc(vector<double>& ta,vector<Field>& evec,
			  int& Nsbt)const{ 

  using namespace FieldExpression;
  const size_t fsize = evec[0].size();

  int Np = Nm_-Nk_;
  CCIO::cout << " -- Nk = " << Nk_<< " Np = "<< Np << endl;
  CCIO::cout << " -- Nm = " << Nm_<< endl;
  CCIO::cout << " -- size of lmd   = " << ta.size() << endl;
  CCIO::cout << " -- size of evec  = " << evec.size() << endl;
  
  assert(Nm_< evec.size() && Nm_< ta.size());
  
  vector<double> tb(Nm_),ta2(Nm_),tb2(Nm_);
  vector<double> Qt(Nm_*Nm_);
  vector<Field>  B(Nm_);
  for(int k=0; k<Nm_; ++k) B[k].resize(fsize);
  
  Field f(fsize);

  vector<int> i_conv;  
  double beta_k;
  int Nconv = -1;

  evec[0] = 1.0/sqrt(fsize);   // initial vector (uniform)

  // Initial Nk steps
  for(int k=0; k<Nk_; ++k) decomp_Lanczos(ta,tb,f,evec,k);

  // Restarting loop begins
  for(int iter=0; iter<Niter_; ++iter){
    CCIO::cout<<"\n iteration = "<< iter << endl;

    for(int k=Nk_; k<Nm_; ++k) decomp_Lanczos(ta,tb,f,evec,k);

    // getting eigenvalues
    for(int k=0; k<Nm_; ++k) ta2[k] = ta[k];
    tb2 = tb;

    diagonalize(ta2,tb2,Qt,Nm_);
    sort_->push(ta2,Nm_);

    // Implicitly shifted QR transformations
    setUnit(Qt);
    for(int p=Nk_; p<Nm_; ++p) qr_decomp(ta,tb,Qt,Nm_,ta2[p],0,Nm_-1);
    
    for(int j=0; j<Nk_+1; ++j){
      B[j] = 0.0;
      for(int k=0; k<Nm_; ++k) B[j] += Qt[j*Nm_+k]*evec[k];
    }
    for(int j=0; j<Nk_+1; ++j) evec[j] = B[j];

    // Compressed vector f and beta(k2)
    f *= Qt[Nm_*(Nk_-1)+Nm_-1];
    f += tb[Nk_-1]*evec[Nk_];
    beta_k = sqrt(f*f);
    CCIO::cout<<" beta(k) = "<<beta_k<<endl;

    evec[Nk_] = (1.0/beta_k)*f;
    tb[Nk_-1] = beta_k;

    // Convergence test
    for(int k=0; k<Nm_; ++k) ta2[k] = ta[k];
    tb2 = tb;

    diagonalize(ta2,tb2,Qt,Nk_);

    for(int j=0; j<Nk_; ++j){
      B[j]=0.0;
      for(int k=0; k<Nk_; ++k) B[j] += Qt[j*Nm_+k]*evec[k];
    }

    i_conv.clear();
    int Kthrs = 0;

    CCIO::cout << setiosflags(ios_base::scientific);
    for(int i=0; i<Nk_; ++i){
      Field v = opr_->mult(B[i]);

      double vnum = B[i]*v;
      double vden = B[i]*B[i];
      ta2[i] = vnum/vden;
      v -= ta2[i]*B[i];
      double vv = v*v;

      CCIO::cout <<" [" <<setw( 3)<< setiosflags(ios_base::right)<<i<<"] ";
      CCIO::cout <<       setw(25)<< setiosflags(ios_base::left) <<ta2[i];
      CCIO::cout <<"  "<< setw(25)<< setiosflags(ios_base::right)<<vv<< endl;

      if(vv<enorm_){
        i_conv.push_back(i);
        if(sort_->saturated(ta2[i],vthrs_)) ++Kthrs;
      }
    }  // i-loop end
    CCIO::cout << resetiosflags(ios_base::scientific);

    CCIO::cout<<" #modes converged: "<<i_conv.size()<<endl;
    if(Kthrs > 0){// (there is a converged eigenvalue larger than Vthrs.)
      Nsbt = i_conv.size()-Kthrs;
      Nconv = iter; 
      break;
    }
    if(iter==Niter_) {
      CCIO::cout<<"\n NOT converged.\n";
      abort();
    }
  } // end of iter loop

  // post process after the conversion
  ta.clear();
  evec.clear();
  for(int i=0; i<i_conv.size(); ++i){
    ta.push_back(ta2[i_conv[i]]);
    evec.push_back(B[i_conv[i]]);
  }
  sort_->push(ta,evec,i_conv.size());

  CCIO::cout << "\n Converged\n Summary :\n";
  CCIO::cout << " -- Iterations  = "<< Nconv  <<"\n";
  CCIO::cout << " -- beta(k)     = "<< beta_k <<"\n";
  CCIO::cout << " -- Kdis        = "<< i_conv.size()<<"\n";
  CCIO::cout << " -- Nsbt        = "<< Nsbt   <<"\n";
}

void EigenModes_IRL::decomp_Lanczos(vector<double>& ta,vector<double>& tb,
				    Field& w,vector<Field>& v,int k)const{
  using namespace FieldExpression;

  w = opr_->mult(v[k]);
  if(k) w -= tb[k-1]*v[k-1];
  
  double ab = v[k]*w;
  w -= ab*v[k];
  ta[k] = ab;

  ab = w*w;
  ab = sqrt(ab);
  w *= 1.0/ab;
  tb[k] = ab;

  if(k) orthogonalize(w,v,k);
  if(k<Nm_-1) v[k+1] = w;
  if(k==Nm_-1) w *= tb[Nm_-1];
}

void EigenModes_IRL::
orthogonalize(Field& w,const vector<Field>& v,int N)const{
  // Schmidt orthogonalization                                   
                                                
  size_t size = w.size();
  assert(size%2 ==0);

  std::slice re(0,size/2,2);
  std::slice im(1,size/2,2);

  for(int j=0; j<N; ++j){
    double prdr = v[j]*w;
    double prdi = v[j].im_prod(w);

    valarray<double> evr(v[j][re]);
    valarray<double> evi(v[j][im]);

    w.add(re, -prdr*evr +prdi*evi);
    w.add(im, -prdr*evi -prdi*evr);
  }
}

void EigenModes_IRL::setUnit(vector<double>& Qt)const{
  int N = sqrt(Qt.size());
  for(int i=0; i<Qt.size(); ++i) Qt[i] = 0.0;
  for(int k=0; k<N; ++k) Qt[k*N+k] = 1.0;
}

void EigenModes_IRL::diagonalize(vector<double>& ta,vector<double>& tb,
				 vector<double>& Qt,int Nk)const{
  setUnit(Qt);
  int Niter = 100*Nm_;
  int kmin = 0;
  int kmax = Nk-1;

  for(int iter=0; iter<Niter; ++iter){
    // determination of 2x2 leading submatrix
    double dsub = ta[kmax]-ta[kmax-1];
    double dd = sqrt(dsub*dsub + 4.0*tb[kmax-1]*tb[kmax-1]);
    double sft = 0.5*(ta[kmax-1]+ta[kmax] +dd*(dsub/fabs(dsub)));

    // transformation
    qr_decomp(ta,tb,Qt,Nk,sft,kmin,kmax);

    // Convergence criterion (redef of kmin and kmax)
    for(int j=kmax; j>kmin; --j){
      double dds = fabs(ta[j-1])+fabs(ta[j]);
      if(fabs(tb[j-1])+dds > dds){
        kmax = j;
	break;
      }
    }
    if(kmax == kmin+1){
      Niter = iter;
      return;
    }
  }
  CCIO::cout << "[QR method] Error - Too many iteration: "<<Niter<<"\n";
  abort();
}

void EigenModes_IRL::qr_decomp(vector<double>& ta,vector<double>& tb,
                                vector<double>& Qt,int Nk,
				double sft,int k_min,int k_max)const{
  // initial process (k = k_min)
  // Givens rotation
  double Fden = 1.0/sqrt((ta[k_min]-sft)*(ta[k_min]-sft) +tb[k_min]*tb[k_min]); 
  double c = (ta[k_min]-sft)*Fden;   //  cos_x
  double s = tb[k_min]*Fden;         //  sin_x

  for(int i=0; i<Nk; ++i){           // accumulation of G(k+1,k)^T
    double Qt_k  = Qt[    k_min*Nm_+i];
    double Qt_kp = Qt[(k_min+1)*Nm_+i];
    Qt[    k_min*Nm_+i] =  c*Qt_k +s*Qt_kp;
    Qt[(k_min+1)*Nm_+i] = -s*Qt_k +c*Qt_kp;
  }

  double a_k = ta[k_min]; 
  double b_k = tb[k_min]; 
  
  ta[k_min] = c*c*a_k +s*s*ta[k_min+1] +2.0*c*s*b_k;
  tb[k_min] = (c*c-s*s)*b_k +c*s*(ta[k_min+1] -a_k);

  double x = s*tb[k_min+1];
    
  ta[k_min+1] = s*s*a_k +c*c*ta[k_min+1] -2.0*c*s*b_k;
  tb[k_min+1] *= c;

  // bulk process 
  for(int k=k_min+1; k<k_max; ++k){
    Fden = 1.0/sqrt(tb[k-1]*tb[k-1] +x*x); // Givens rotation
    c = tb[k-1]*Fden;                      //  cos_x
    s = x*Fden;                            //  sin_x

    for(int i=0; i<Nk; ++i){           // accumulation of G(k+1,k)^T
      double Qt_k  = Qt[    k*Nm_+i];
      double Qt_kp = Qt[(k+1)*Nm_+i];
      Qt[    k*Nm_+i] =  c*Qt_k +s*Qt_kp;
      Qt[(k+1)*Nm_+i] = -s*Qt_k +c*Qt_kp;
    }
    a_k = ta[k]; 
    b_k = tb[k]; 

    // G(k+1,k)^T A G(k+1,k)
    tb[k-1]= c*tb[k-1] +s*x;                    // k-1,k
    ta[k  ]= c*c*a_k +s*s*ta[k+1] +2.0*c*s*b_k; // k,k
    tb[k  ]= (c*c-s*s)*b_k +c*s*(ta[k+1] -a_k); // k,k+1

    x = s*tb[k+1];
    
    ta[k+1]= s*s*a_k +c*c*ta[k+1] -2.0*c*s*b_k; // k+1,k+1
    tb[k+1]*= c;                                // k+1,k+2
  }
}

