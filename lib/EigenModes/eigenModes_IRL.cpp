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
#include "commonPrms.h"
using namespace std;

void EigenModes_IRL::calc(vector<double>& ta,vector<Field>& V,int& Nin)const{ 

  using namespace FieldExpression;
  const size_t fsize = opr_->fsize();

  Field f(fsize);

  int Np = Nm_-Nk_;  
  assert(Np>0);
  CCIO::cout <<"Nk = "<< Nk_<<" Np = "<< Np<<" Nm = "<<Nm_<< endl;
  
  ta.resize(Nm_);  V.resize(Nm_);
  vector<double> tb(Nm_);

  V[0] = 1.0/sqrt(double(fsize));   // initial vector (uniform)
  lanczos_init(ta,tb,V);            // initial Lanczos-decomp
  //for(int k=0; k<Nk_; ++k) step(Nm_,k,ta,tb,V,f);

  CCIO::cout<<"V[0][0]="<<V[0][0]<<endl;

  // Implicit Restarting
  vector<double> tta(Nm_),ttb(Nm_);
  vector<double> Qt(Nm_*Nm_);

  vector<Field>  Vp(Nk_+1);
  for(int k=0; k<Nk_+1; ++k) Vp[k].resize(fsize);

  CCIO::cout<<"lanczos decomp(init)"<<endl;
  for(int j=0;j<ta.size();++j) CCIO::cout<<j<<" "<<ta[j]<<" "<<tb[j]<<endl;

  vector<int> i_conv;  
  int Nconv = -1;

  for(int iter=0; iter<Niter_; ++iter){
    CCIO::cout<<"\n iteration "<< iter << endl;

    //    for(int k=Nk_; k<Nm_; ++k) step(Nm_,k,ta,tb,V,f);
    //    f *= tb[Nm_-1];

    lanczos_ext(ta,tb,V,f);         // extended Lanczos-decomp
    tta = ta;  ttb = tb;            // getting shifts 

    CCIO::cout<<"lanczos decomp(extd)"<<endl;
    for(int j=0;j<ta.size();++j) CCIO::cout<<j<<" "<<ta[j]<<" "<<tb[j]<<endl;

    diagonalize(tta,ttb,Qt,Nm_); 
    /*
    CCIO::cout<<"before sort"<<endl;
    for(int j=0;j<tta.size();++j) CCIO::cout<<j<<" "<<tta[j]<<endl;
    */
    sort_->push(tta,Nm_);          // sort by the absolute values
    /*
    CCIO::cout<<"after sort"<<endl;
    for(int j=0;j<tta.size();++j) CCIO::cout<<j<<" "<<tta[j]<<endl;
    */
    // QR transformations with implicit shifts tta[p]
    setUnit(Qt);
    for(int p=Nk_; p<Nm_; ++p) {
      QRdecomp_Givens(ta,tb,Qt,Nm_,tta[p],0,Nm_-1);
      CCIO::cout<<"p="<<p<<endl;
      for(int j=0;j<ta.size();++j) CCIO::cout<<j<<" "<<ta[j]<<" "<<tb[j]<<endl;
    }

    for(int i=0; i<Nk_+1; ++i){         // getting V+
      Vp[i] = 0.0;
      for(int j=0; j<Nm_; ++j) Vp[i] += Qt[i*Nm_+j]*V[j];
    }
    for(int i=0; i<Nk_+1; ++i) V[i] = Vp[i];

    f *= Qt[Nm_*(Nk_-1)+Nm_-1];        // getting f+
    f += tb[Nk_-1]*V[Nk_];
    //tb[Nk_-1] = sqrt(f*f);
    //V[Nk_] = 1.0/tb[Nk_-1]*f;

    // Convergence test
    tta = ta;  ttb = tb;
    diagonalize(tta,ttb,Qt,Nk_); // Qt contains Nk_ eigenvectors

    i_conv.clear();
    int Nover = 0; // Num of converged eigenvalues beyond Vthrs.

    CCIO::cout << setiosflags(ios_base::scientific);

    for(int i=0; i<Nk_; ++i){
      Vp[i] = 0.0;                    // eigenvectors
      for(int j=0; j<Nk_; ++j) Vp[i] += Qt[i*Nm_+j]*V[j];  

      Field Av = opr_->mult(Vp[i]);

      double vAv = Vp[i]*Av;        
      double vv = Vp[i]*Vp[i];          
      tta[i] = vAv/vv;                // eigenvalues
      Av -= tta[i]*Vp[i];
      double res = Av*Av;

      CCIO::cout<<" ["<<setw( 3)<<setiosflags(ios_base::right)<<i<<"] ";
      CCIO::cout<<      setw(25)<<setiosflags(ios_base::left) <<tta[i];
      CCIO::cout<<"  "<<setw(25)<<setiosflags(ios_base::right)<<res<<endl;

      if(res<enorm_){
        i_conv.push_back(i);
        if(sort_->beyond_thrs(tta[i],vthrs_)) ++Nover;
      }
    }  // i-loop end
    CCIO::cout << resetiosflags(ios_base::scientific);

    CCIO::cout<<" #converged modes= "<<i_conv.size()<<endl;
    if(Nover > 0){ 
      Nin = i_conv.size()-Nover; 
      Nconv = iter; 
      break;
    }
    if(iter==Niter_) {
      CCIO::cout<<"\n NOT converged.\n";
      abort();
    }
  } // end of iter loop

  // post process after the conversion
  ta.clear(); V.clear();
  
  for(int i=0; i<i_conv.size(); ++i){
    ta.push_back(tta[i_conv[i]]);
    V.push_back(Vp[i_conv[i]]);
  }

  sort_->push(ta,V,i_conv.size());
  double beta_pk = sqrt(f*f);
  CCIO::cout << "\n Converged\n Summary :\n";
  CCIO::cout << " -- Iterations  = "<< Nconv        <<"\n";
  //CCIO::cout << " -- beta+_k     = "<< tb[Nk_-1]    <<"\n";
  CCIO::cout << " -- beta+_k     = "<< beta_pk    <<"\n";
  CCIO::cout << " -- N_eigen     = "<< i_conv.size()<<"\n";
  CCIO::cout << " -- N_inner     = "<< Nin          <<"\n";
}

void EigenModes_IRL::
lanczos_init(vector<double>& ta,vector<double>& tb,vector<Field>& V)const{
  using namespace FieldExpression;
  
  Field f = opr_->mult(V[0]);
  double ab = V[0]*f;
  f -= ab*V[0];
  ta[0] = ab;

  ab = f*f;
  tb[0] = sqrt(ab);
  V[1] = 1.0/tb[0]*f;
  
  for(int k=1; k<Nk_;++k){
    f = opr_->mult(V[k]);
    f -= tb[k-1]*V[k-1];
  
    ab = V[k]*f;
    f -= ab*V[k];
    ta[k] = ab;

    ab = f*f;
    tb[k] = sqrt(ab);
    f *= 1.0/tb[k]; 

    orthogonalize(f,V,k);
    ta[k] += V[k]*f;
    tb[k-1] += V[k-1]*f;

    V[k+1] = f;
  }
}

void EigenModes_IRL::lanczos_ext(vector<double>& ta,vector<double>& tb,
				 vector<Field>& V,Field& f)const{
  using namespace FieldExpression;
  
  for(int k=Nk_; k<Nm_; ++k){
    f = opr_->mult(V[k]);
    f -= tb[k-1]*V[k-1];
  
    double ab = V[k]*f;
    f -= ab*V[k];
    ta[k] = ab;

    ab = f*f;
    tb[k] = sqrt(ab);
    f *= 1.0/tb[k];

    orthogonalize(f,V,k);
    
    if(k!= Nm_-1) V[k+1] = f;
  }
  f *= tb[Nm_-1];
}


void EigenModes_IRL::
orthogonalize(Field& f,const vector<Field>& V,int k)const{
  // Schmidt orthogonalization                                   
                                                
  size_t size = f.size();
  assert(size%2 ==0);

  std::slice re(0,size/2,2);
  std::slice im(1,size/2,2);

  vector<double> sr(size/2);
  vector<double> si(size/2);

  for(int j=0; j<k; ++j){
    sr[j] = V[j]*f;
    si[j] = V[j].im_prod(f);
  }
  for(int j=0; j<k; ++j){
    f.add(re, -sr[j]*V[j][re] +si[j]*V[j][im]);
    f.add(im, -sr[j]*V[j][im] -si[j]*V[j][re]);
  }
}
/*
void EigenModes_IRL::
orthogonalize(Field& f,const vector<Field>& V,int k)const{
  // Schmidt orthogonalization                                   
                                                
  size_t size = f.size();
  assert(size%2 ==0);

  std::slice re(0,size/2,2);
  std::slice im(1,size/2,2);

  for(int j=0; j<k; ++j){
    double prdr = V[j]*f;
    double prdi = V[j].im_prod(f);

//    valarray<double> fre = f[re];
//    valarray<double> fim = f[im];
    
//    fre -= prdr*V[j][re]; fre += prdi*V[j][im];
//    fim -= prdr*V[j][im]; fim -= prdi*V[j][re];

//    f.set(re,fre);
//    f.set(im,fim);

    valarray<double> vr(V[j][re]);
    valarray<double> vi(V[j][im]);

    f.add(re, -prdr*vr +prdi*vi);
    f.add(im, -prdr*vi -prdi*vr);
  }
}
*/

void EigenModes_IRL::setUnit(vector<double>& Qt)const{
  for(int i=0; i<Nm_*Nm_; ++i) Qt[i] = 0.0;
  for(int k=0; k<Nm_; ++k) Qt[k*Nm_+k] = 1.0;
}

void EigenModes_IRL::diagonalize(vector<double>& ta,vector<double>& tb,
				 vector<double>& Qt,int Nk)const{
  setUnit(Qt);
  int Niter = 100*Nm_;
  int kmin = 0;
  int kmax = Nk-1;

  for(int iter=0; iter<Niter; ++iter){
    // looking at the rightmost 2x2 submatrix
    double sb = ta[kmax]-ta[kmax-1];
    double dd = sqrt(sb*sb + 4.0*tb[kmax-1]*tb[kmax-1]);
    double sft = 0.5*(ta[kmax-1]+ta[kmax] +dd*(sb/fabs(sb)));

    // transformation
    QRdecomp_Givens(ta,tb,Qt,Nk,sft,kmin,kmax);
    //for(int j=0; j<ta.size(); ++j) CCIO::cout<<j<<" "<<ta[j]<<endl;

    // Convergence criterion (redef of kmax)
    if(kmax == kmin+1){
      Niter = iter;
      return;
    }
    for(int j=kmax; j>kmin; --j){
      double dds = fabs(ta[j-1])+fabs(ta[j]);
      if(fabs(tb[j-1])+dds > dds){
        kmax = j;
	break;
      }
    }
  }
  CCIO::cout << "[QR method] Error - Too many iteration: "<<Niter<<"\n";
  abort();
}

void EigenModes_IRL::QRdecomp_Givens(vector<double>& ta,vector<double>& tb,
				     vector<double>& Qt,int Nk,
				     double sft,int k_min,int k_max)const{
  // initial process (k = k_min)
  // Givens rotation
  double dn = 1.0/sqrt((ta[k_min]-sft)*(ta[k_min]-sft)+tb[k_min]*tb[k_min]); 
  double c = (ta[k_min]-sft)*dn;   //  cos_x
  double s = tb[k_min]*dn;         //  sin_x

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
    dn = 1.0/sqrt(tb[k-1]*tb[k-1] +x*x); // Givens rotation
    c = tb[k-1]*dn;                      //  cos_x
    s = x*dn;                            //  sin_x

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


//////////////////////////////////////////////////////////////////////////
/*
void EigenModes_IRL::step(int Nm, int k, vector<double>& TDa,
			  vector<double>& TDb,
			  vector<Field>& vk, Field& w)const{
  using namespace FieldExpression;
  if(k>=Nm){
    // pprintf("k larger than Nm: irlegal.\n");
    abort();
  }else if(k==0){  // Initial step

    w = opr_->mult(vk[k]);
    double alph = vk[k] * w;

    w -= alph * vk[k];
    double beta = w * w;
    beta = sqrt(beta);
    double betar = 1.0/beta;
    vk[k+1] = betar * w;
    TDa[k] = alph;
    TDb[k] = beta;

  }else{   // Iteration step

    w = opr_->mult(vk[k]);
    w -= TDb[k-1] * vk[k-1];
    double alph = vk[k] * w;

    w -= alph * vk[k];
    double beta = w * w;
    beta = sqrt(beta);
    double betar = 1.0/beta;
    w *= betar;
    TDa[k] = alph;
    TDb[k] = beta;

    orthogonalize(w, vk, k);

    if(k < Nm-1) vk[k+1] = w;

  };

};
*/
/*
void EigenModes_IRL::
lanczos_init(vector<double>& TDa,
	     vector<double>& TDb,
	     vector<Field>& vk) const{

  using namespace FieldExpression;
  Field w = opr_->mult(vk[0]);
  double alph = vk[0] * w;
    
  w -= alph * vk[0];
  double beta = w * w;
  beta = sqrt(beta);
  double betar = 1.0/beta;
  vk[1] = betar * w;
  TDa[0] = alph;
  TDb[0] = beta;

  for(int k=1;k<Nk_;++k){

    w = opr_->mult(vk[k]);
    w -= TDb[k-1] * vk[k-1];
    double alph = vk[k] * w;

    w -= alph * vk[k];
    double beta = w * w;
    beta = sqrt(beta);
    double betar = 1.0/beta;
    w *= betar;
    TDa[k] = alph;
    TDb[k] = beta;

    orthogonalize(w, vk, k);

    vk[k+1] = w;

  };

};
*/

/*
void EigenModes_IRL::lanczos_ext(vector<double>& TDa,
				 vector<double>& TDb,
				 vector<Field>& vk, Field& w)const{
  using namespace FieldExpression;
  for(int k=Nk_;k<Nm_;++k){
    w = opr_->mult(vk[k]);
    w -= TDb[k-1] * vk[k-1];
    double alph = vk[k] * w;

    w -= alph * vk[k];
    double beta = w * w;
    beta = sqrt(beta);
    double betar = 1.0/beta;
    w *= betar;
    TDa[k] = alph;
    TDb[k] = beta;
    
    orthogonalize(w, vk, k);

    if(k < Nm_-1) vk[k+1] = w;
  }
  w *=TDb[Nm_-1];
};
*/

/*
void EigenModes_IRL::
orthogonalize(Field& w, const vector<Field>& vk, int k)const{
  // Schmidt orthogonalization


//  for(int j=0; j<k; ++j){
//    double prd = vk[j] * w;
//    w -= prd * vk[j];
//    // Originally this operation is in complex.
//    // It is NOT sufficient with real operation, and
//    // complex version is needed.
//  };
//

  //if( (fmt.Nin()%2)!=0 ) abort();
  Format::Format_F fmt(CommonPrms::instance()->Nvol());
  int Nin  = fmt.Nin();
  int Nvol = fmt.Nvol();
  int Nex  = fmt.Nex();
  double prdr, prdi;
  double vr, vi;

  
  for(int j=0; j<k; ++j){

    prdr = 0.0;
    prdi = 0.0;
    for(int ex=0; ex<Nex; ++ex){
     for(int iv=0; iv<Nvol; ++iv){
      for(int in=0; in<Nin; in+=2){
        prdr += vk[j][in+Nin*(iv+Nvol*ex)]*w[in+Nin*(iv+Nvol*ex)]
	       +vk[j][in+1+Nin*(iv+Nvol*ex)]*w[in+1+Nin*(iv+Nvol*ex)];
        prdi += vk[j][in+Nin*(iv+Nvol*ex)]*w[in+1+Nin*(iv+Nvol*ex)]
	       -vk[j][in+1+Nin*(iv+Nvol*ex)]*w[in+Nin*(iv+Nvol*ex)];
      };
     };
    };
    prdr = Communicator::instance()->reduce_sum(prdr);
    prdi = Communicator::instance()->reduce_sum(prdi);

    for(int ex=0; ex<Nex; ++ex){
     for(int iv=0; iv<Nvol; ++iv){
      for(int in=0; in<Nin; in+=2){
        vr = w[in+Nin*(iv+Nvol*ex)] -prdr*vk[j][in+Nin*(iv+Nvol*ex)]
                          	  + prdi*vk[j][in+1+Nin*(iv+Nvol*ex)];
        vi = w[in+1+Nin*(iv+Nvol*ex)]- prdr*vk[j][in+1+Nin*(iv+Nvol*ex)]
                               - prdi*vk[j][in+Nin*(iv+Nvol*ex)];
        w.set(in+Nin*(iv+Nvol*ex),vr);
        w.set(in+1+Nin*(iv+Nvol*ex),vi);
      };
     };
    };

  };

};

*/


