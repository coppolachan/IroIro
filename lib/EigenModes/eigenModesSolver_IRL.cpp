/*! @file eigenModesSolver_IRL.cpp
 *  @brief implementations of EigenModesSolver_IRL class  
 */
#include "eigenModesSolver_IRL.hpp"
#include "eigenSorter.hpp"
#include "Fopr/fopr.h"
#include "Fields/field_expressions.hpp"
#include "include/messages_macros.hpp"
#include <cassert>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

void EigenModesSolver_IRL::
calc(vector<double>& ta,vector<Field>& V,int& Neigen)const{ 

  _Message(DEBUG_VERB_LEVEL, "EigenModesSolver_IRL::calc\n");
  using namespace FieldExpression;
  const size_t fsize = opr_->fsize();

  Field f(fsize);

  int Np = Nm_-Nk_;  
  CCIO::cout <<"Nk = "<< Nk_<<" Np = "<< Np<<" Nm = "<<Nm_<<"\n";

  ta.resize(Nm_); 
  V.assign(Nm_,Field(fsize,1.0));

  double nv = V[0].norm();
  V[0] /= nv;                       /*!< @brief initial vector (uniform)*/

  vector<double> tb(Nm_);
  lanczos_init(ta,tb,V);            /*!< @brief initial Lanczos-decomp */

  vector<double> tta(Nm_),ttb(Nm_);
  vector<double> Qt(Nm_*Nm_);

  vector<Field>  Vp(Nk_+1);
  for(int k=0; k<Nk_+1; ++k) Vp[k].resize(fsize);

  vector<int> idx(Nk_);  
  int Iconv = -1; /*!<@brief Num of iterations until convergence */
  int Ncert = 0;  /*!<@brief Num of the certified eigenmodes */

  for(int iter=0; iter<Niter_; ++iter){
    CCIO::cout<<"\n iteration "<< iter << endl;

    /****** Restarting procedures ******/
    lanczos_ext(ta,tb,V,f);        /*!< @brief extended Lanczos-decomp*/
    tta = ta;  ttb = tb;           /*!< @brief getting shifts */

    diagonalize(tta,ttb,Qt,Nm_); 
    esorter_->push(tta,Nm_);          /*!< @brief sort by the absolute values*/

    /// QR transformations with implicit shifts tta[p]
    setUnit(Qt);
    for(int p=Nk_; p<Nm_; ++p) QRfact_Givens(ta,tb,Qt,Nm_,tta[p],0,Nm_-1);

    for(int i=0; i<Nk_+1; ++i){         // getting V+
      Vp[i] = 0.0;
      for(int j=0; j<Nm_; ++j) Vp[i] += Qt[i*Nm_+j]*V[j];
    }
    for(int i=0; i<Nk_+1; ++i) V[i] = Vp[i];

    f *= Qt[Nm_*(Nk_-1)+Nm_-1];        // getting f+
    f += tb[Nk_-1]*V[Nk_];

    tb[Nk_-1] = sqrt(f*f);
    V[Nk_] = 1.0/tb[Nk_-1]*f;

    /******  Convergence test  ******/
    tta = ta;  ttb = tb;
    diagonalize(tta,ttb,Qt,Nk_); /*!< @brief Qt contains Nk_ eigenvectors */
    
    vector<double> res(Nm_);  
    int Nover = 0; /*!< @brief Num of converged eigenvalues beyond thrs.*/

    CCIO::cout << setiosflags(ios_base::scientific);

    for(int i=0; i<Nk_; ++i){
      idx[i] = i;                     /*!< @brief eigen-indices */

      Vp[i] = 0.0;                    /*!< @brief eigenvectors */
      for(int j=0; j<Nk_; ++j) Vp[i] += Qt[i*Nm_+j]*V[j];  

      Field Av = opr_->mult(Vp[i]);
      tta[i] = Vp[i]*Av;
      Av -= tta[i]*Vp[i];
      res[i] = Av.norm();

      CCIO::cout<<" ["<<setw( 3)<<setiosflags(ios_base::right)<<i<<"] ";
      CCIO::cout<<      setw(25)<<setiosflags(ios_base::left) <<tta[i];
      CCIO::cout<<"  "<<setw(25)<<setiosflags(ios_base::right)<<res[i]<<endl;
    }

    Ncert =0;
    esorter_->push(tta,idx,Nk_);

    for(int i=0; i<Nk_; ++i){
      if(res[i]<prec_){  /*!<@brief counting converged eigenmodes */
	Ncert++;
        if(esorter_->beyond_thrs(tta[i],opr_)) ++Nover; 
      }else{
	break;
      }
    }  
    CCIO::cout << resetiosflags(ios_base::scientific);
    CCIO::cout<<" #converged modes= "<<Ncert <<endl;
     
    /*!@brief  condition of termination: the eigenvalue exceeds the threshold*/
    if(Nover > 0){         
      CCIO::cout<<"All eigenmodes upper/under the threshold are obtained.\n";
      Iconv = iter; 
      Neigen = Ncert -Nover; 
      break;
    }
    /*!@brief condition of termination: #eigenmods exceeds the threshold number*/
    if(Ncert >= Nthrs_){ 
      CCIO::cout<<"Desired number of eigenmodes are obtained.\n";
      Iconv = iter;
      Neigen = Nthrs_; 
      break;
    }
    
    if(iter==Niter_-1){
      CCIO::cout<<"Reached the max iteration count.\n";
      Neigen = -Ncert;
      break;
    } 
  }// end of iter loop

  
  /*** post process after the conversion ***/
  ta.clear();
  V.clear();

  
  for(int i=0; i < Neigen; i++){
    ta.push_back(tta[i]);
    V.push_back(Vp[idx[i]]);
  }

  if(Neigen > 0){
    CCIO::cout << "\n Converged\n Summary :\n";
    CCIO::cout << " -- Iterations  = "<< Iconv     <<"\n";
    CCIO::cout << " -- beta+_k     = "<< tb[Nk_-1] <<"\n";
    CCIO::cout << " -- #converged  = "<< Ncert     <<"\n";
    CCIO::cout << " -- #eigenmodes = "<< Neigen    <<"\n";
  }
}

void EigenModesSolver_IRL::
lanczos_init(vector<double>& ta,vector<double>& tb,vector<Field>& V)const{
  _Message(DEBUG_VERB_LEVEL, "EigenModesSolver_IRL::lanczos_init\n");
  using namespace FieldExpression;

  Field f = opr_->mult(V[0]);

  double ab = V[0]*f;
  f -= ab*V[0];
  ta[0] = ab;
  ab = f*f;
  tb[0] = sqrt(ab);
  V[1] = 1.0/tb[0]*f;


  _Message(DEBUG_VERB_LEVEL, "EigenModesSolver_IRL::lanczos_init Starting loop\n");
  for(int k=1; k<Nk_;++k){
    _Message(DEBUG_VERB_LEVEL, "EigenModesSolver_IRL::lanczos_init loop k="<<k << " mult\n" );
    f = opr_->mult(V[k]);
    _Message(DEBUG_VERB_LEVEL, "EigenModesSolver_IRL::lanczos_init loop k="<<k << " after mult\n" );
    f -= tb[k-1]*V[k-1];

    ab = V[k]*f;
    f -= ab*V[k];
    ta[k] = ab;

    ab = f*f;
    tb[k] = sqrt(ab);
    f *= 1.0/tb[k]; 

    orthogonalize(f,V,k); /*!< classical Gram-Schmidt orthogonalization */
    ta[k] += V[k]*f;      /*!<@brief DGKS correction*/
    tb[k-1] += V[k-1]*f;

    V[k+1] = f;
  }
}

void EigenModesSolver_IRL::lanczos_ext(vector<double>& ta,vector<double>& tb,
				       vector<Field>& V,Field& f)const{
  using namespace FieldExpression;
  _Message(DEBUG_VERB_LEVEL, "EigenModesSolver_IRL::lanczos_ext\n");  

  for(int k=Nk_; k<Nm_; ++k){
    f = opr_->mult(V[k]);
    f -= tb[k-1]*V[k-1];
  
    double ab = V[k]*f;
    f -= ab*V[k];
    ta[k] = ab;

    ab = f*f;
    tb[k] = sqrt(ab);
    f *= 1.0/tb[k];

    orthogonalize(f,V,k); /*!<@brief classical Gram-Schmidt orthogonalization*/
    ta[k] += V[k]*f;      /*!<@brief DGKS correlction*/
    tb[k-1] += V[k-1]*f;

    if(k!= Nm_-1) V[k+1] = f;
  }
  f *= tb[Nm_-1];
}

// classical Gram-Schmidt orthogonalization
void EigenModesSolver_IRL::
orthogonalize(Field& f,const vector<Field>& V,int k)const{
  _Message(DEBUG_VERB_LEVEL, "EigenModesSolver_IRL::orthogonalize\n");  

  size_t size = f.size();
  assert(size%2 ==0);

  std::slice re(0,size/2,2);
  std::slice im(1,size/2,2);

  vector<double> sr(k);
  vector<double> si(k);

  for(int j=0; j<k; ++j){
    sr[j] = V[j]*f;
    si[j] = V[j].im_prod(f);
  }
  for(int j=0; j<k; ++j){
    f.add(re, -sr[j]*V[j][re] +si[j]*V[j][im]);
    f.add(im, -sr[j]*V[j][im] -si[j]*V[j][re]);
  }
}

void EigenModesSolver_IRL::setUnit(vector<double>& Qt)const{
  for(int i=0; i<Nm_*Nm_; ++i) Qt[i] = 0.0;
  for(int k=0; k<Nm_; ++k) Qt[k*Nm_+k] = 1.0;
}

void EigenModesSolver_IRL::diagonalize(vector<double>& ta,vector<double>& tb,
				       vector<double>& Qt,int Nk)const{

  _Message(DEBUG_VERB_LEVEL, "EigenModesSolver_IRL::diagonalize\n");  
  setUnit(Qt);
  int Niter = 100*Nm_;
  int kmin = 0, kmax = Nk-1;

  for(int iter=0; iter<Niter; ++iter){

    double sft = 0.5*(ta[kmax]+ta[kmax-1]);
    double dif = 0.5*(ta[kmax]-ta[kmax-1]);
    sft += dif/fabs(dif)*sqrt(dif*dif +tb[kmax-1]*tb[kmax-1]);
    
    QRfact_Givens(ta,tb,Qt,Nk,sft,kmin,kmax); // transformation

    for(int j=kmax; j>kmin; --j){   
      double dds = fabs(ta[j-1])+fabs(ta[j]);
      if(fabs(tb[j-1])+dds > dds){
        kmax = j;
	break;
      }
      if(j==kmin+1){/*!< Convergence criterion */
	Niter = iter;
	return; 
      }
    }
  }
  CCIO::cout<<"[QR method] reached to maximum iterations: "<<Niter<<"\n";
  CCIO::cout<<"failed at kmax="<<kmax<<"\n";
  abort();
}

void EigenModesSolver_IRL::QRfact_Givens(vector<double>& ta,vector<double>& tb,
					 vector<double>& Qt,int Nk,
					 double sft,int k_min,int k_max)const{
  _Message(DEBUG_VERB_LEVEL, "EigenModesSolver_IRL::QRfact_Givens\n");  
  // k = k_min
  double r = tb[k_min]/(ta[k_min]-sft);
  double c = 1.0/sqrt(1.0+r*r);   //  cos_x
  double s =   r/sqrt(1.0+r*r);   //  sin_x
  
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

  // k >k_min
  for(int k=k_min+1; k<k_max; ++k){
    r = x/tb[k-1];           /*!< Givens rotation */
    c = 1.0/sqrt(1.0+r*r);   //  cos_x
    s =   r/sqrt(1.0+r*r);   //  sin_x

    for(int i=0; i<Nk; ++i){           /*!< accumulation of G(k+1,k)^T */
      double Qt_k  = Qt[    k*Nm_+i];
      double Qt_kp = Qt[(k+1)*Nm_+i];
      Qt[    k*Nm_+i] =  c*Qt_k +s*Qt_kp;
      Qt[(k+1)*Nm_+i] = -s*Qt_k +c*Qt_kp;
    }
    a_k = ta[k]; 
    b_k = tb[k]; 

    // unitary transformation G(k+1,k)^T A G(k+1,k)
    tb[k-1]= c*tb[k-1] +s*x;                    // k-1,k
    ta[k  ]= c*c*a_k +s*s*ta[k+1] +2.0*c*s*b_k; // k,k
    tb[k  ]= (c*c-s*s)*b_k +c*s*(ta[k+1] -a_k); // k,k+1

    x = s*tb[k+1];
    
    ta[k+1]= s*s*a_k +c*c*ta[k+1] -2.0*c*s*b_k; // k+1,k+1
    tb[k+1]*= c;                                // k+1,k+2
  }
}

