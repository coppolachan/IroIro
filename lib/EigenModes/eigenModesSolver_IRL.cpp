/*! @file eigenModesSolver_IRL.cpp
 *  @brief implementations of EigenModesSolver_IRL class  
 Time-stamp: <2014-10-09 15:03:53 cossu>
 */
#include "eigenModesSolver_IRL.hpp"
#include "eigenSorter.hpp"
#include "subSpaceProjector.hpp"
#include "Fopr/fopr.h"
#include "Fields/field_expressions.hpp"
#include "include/messages_macros.hpp"
#include <cassert>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

void EigenModesSolver_IRL::calc(vector<double>& ta,vector<Field>& V,int& Neigen,
				vector<Field>* exvec)const{ 
  
  _Message(DEBUG_VERB_LEVEL, "EigenModesSolver_IRL::calc\n");

  using namespace FieldExpression;
  const size_t fsize = opr_->fsize();

  CCIO::cout <<"Nk = "<< Nk_<<" Np = "<< Nm_-Nk_<<" Nm = "<<Nm_<<"\n";
  CCIO::cout<<"Ritz pair with res > "<< ngPrec_*prec_<<" is ignored\n";

  V.assign(Nm_,Field(fsize,1.0));

  double nv = V[0].norm();
  V[0] /= nv;                   /*!< @brief initial vector (uniform)*/

  Field f = opr_->mult(V[0]);

  ta.resize(Nm_); 
  vector<double> tb(Nm_); 

  ta[0] = V[0]*f;
  f -= ta[0]*V[0];
  lanczos_ext(ta,tb,V,f,1,Nk_); /*!< @brief initial Lanczos-decomp */

  vector<double> Qt(Nm_*Nm_);
  vector<double> eval(Nk_);  
  vector<Field> Vp(Nk_,Field(fsize));

  int Iconv = -1;               /*!<@brief Num of iterations until convergence*/
  vector<pair<int,int> > Icert; /*!<@brief index of the certified eigenmodes */

  ///////////////* iteration */////////////////////

  for(int iter=0; iter<Niter_; ++iter){
    CCIO::cout<<"\n iteration "<< iter << endl;

    //// iteration step ////
    if(exvec){                           /*!<@brief exvec_ is projected out */
      for(int k=0; k<Nk_; ++k){  
	SubSpace::projectOut(V[Nk_],V[k],*exvec);
	V[k] = V[Nk_];                   /*!<@brief V[Nk_] is used as buffer*/
      }
    }
    
    lanczos_ext(ta,tb,V,f,Nk_,Nm_);      /*!<@brief Nm-step Lanczos-decomp */

    vector<double> tta = ta;             /*!<@brief tta: eigenvals of H(m) */
    vector<double> ttb = tb;
    int Ndiag = diagonalize(tta,ttb,Qt,Nm_);

    esorter_->push(tta,Nm_);            /*!<@brief sort of absolute values */

    setUnit(Qt);                         /*!<@brief QR-shifts tta[p] */
    for(int p=Nk_; p<Nm_; ++p) QRfact_Givens(ta,tb,Qt,Nm_,tta[p],0,Nm_-1);

    for(int i=0; i<Nk_; ++i){               /*!<@brief getting V+(Nk) */
      Vp[i] = 0.0;
      int np = Nm_-Nk_+i+1;
      for(int j=0; j<np; ++j) Vp[i] += Qt[i*Nm_+j]*V[j];
    }

    f *= Qt[Nm_*(Nk_-1)+Nm_-1];             /*!<@brief getting f+(Nk) */
    for(int j=0; j<Nm_; ++j) f += tb[Nk_-1]*Qt[Nk_*Nm_+j]*V[j];
    for(int i=0; i<Nk_; ++i) V[i] = Vp[i];

    ///  Convergence test ///
    tta = ta;  ttb = tb;          
    Ndiag = diagonalize(tta,ttb,Qt,Nk_);   /*!<@brief Qt*H*Q = tta */

    vector<double> res(Nk_);               /*!<@brief residual     */
    vector<int> idx(Nk_);                  /*!<@brief sorted index */

    for(int i=0; i<Nk_; ++i){
      idx[i] = i;                          /*!<@brief eigen-index  */
      eval[i] = tta[i];                    /*!<@brief Ritz value   */
      Vp[i] = 0.0;                         /*!<@brief Ritz vector  */
      for(int j=0; j<Nk_; ++j) Vp[i] += Qt[i*Nm_+j]*V[j];  
      res[i] = fabs(Qt[i*Nm_+Nk_-1]);
      res[i] *= f.norm();
      //CCIO::cout<<"eval["<<i<<"]="<<eval[i]<<" res="<<res[i]<<"\n";
    }
    esorter_->push(eval,idx,Nk_);

    CCIO::cout << setiosflags(ios_base::scientific);
    for(int i=0; i<Nk_; ++i)
      CCIO::cout<<" ["<<setw( 3)<<setiosflags(ios_base::right)<<i<<"] "
		<<      setw(25)<<setiosflags(ios_base::left) <<eval[i]
		<<"  "<<setw(25)<<setiosflags(ios_base::right)<<res[idx[i]]<<"\n";

    int Nover = 0;              /*!<@brief Num of converged modes beyond thrs.*/
    Icert.clear();

    for(int i=0; i<Nk_; ++i){
      if(res[idx[i]]<prec_){    /*!<@brief counting converged modes */
	Icert.push_back(pair<int,int>(i,idx[i]));
        if(esorter_->beyond_thrs(eval[i],opr_)) ++Nover; 
      }else{
	if(res[idx[i]] >ngPrec_*prec_) continue;
	else                           break;
      }
    }  
    /*
    for(int i=0;i<Icert.size();++i) 
      CCIO::cout<<"i="<<i
		<<" Icert.first="<<Icert[i].first
		<<" Icert.second="<<Icert[i].second<<"\n";
    */
    CCIO::cout << resetiosflags(ios_base::scientific);
    CCIO::cout<<" #converged modes= "<<Icert.size() <<"\n";
     
    if(Nover > 0){           /*!<@brief case1: having eigenvalue beyond thrs. */
      Iconv = iter; 
      Neigen = Icert.size() -Nover; 
      CCIO::cout<<"All eigenmodes upper/under the threshold are obtained.\n";
      break;
    }

    if(Icert.size()>=Nthrs_){/*!<@brief case2: having more modes than Nthrs */
      Iconv = iter;
      Neigen = Nthrs_; 
      CCIO::cout<<"Desired number of eigenmodes are obtained.\n";
      break;
    }

    if(iter==Niter_-1){      /*!<@brief case3: reached to max-iteration */
      Neigen = -Icert.size();
      CCIO::cout<<"Reached the max iteration count.\n";
      break;
    } 
  }// end of iter loop

  /////////* post process after convergence *///////////////
  ta.clear(); V.clear();
  
  for(int i=0; i<Neigen; ++i){
    ta.push_back(eval[Icert[i].first]);
    V.push_back(Vp[Icert[i].second]);
  }

  if(Neigen > 0){
    CCIO::cout << "\n Converged\n Summary :\n";
    CCIO::cout << " -- Iterations  = "<< Iconv       <<"\n";
    CCIO::cout << " -- #converged  = "<< Icert.size()<<"\n";
    CCIO::cout << " -- #eigenmodes = "<< Neigen      <<"\n";
  }
}

void EigenModesSolver_IRL::
lanczos_ext(vector<double>& ta,vector<double>& tb,
	    vector<Field>& V,Field& f,int ini,int fin)const{
  using namespace FieldExpression;

  tb[ini-1] = f.norm();
  V[ini] = 1.0/tb[ini-1]*f;
  
  for(int k=ini; k<fin; ++k){
    f = opr_->mult(V[k]);
    f -= tb[k-1]*V[k-1];
  
    ta[k] = V[k]*f;
    f -= ta[k]*V[k];
    tb[k] = f.norm();
    f /= tb[k];

    orthogonalize(f,V,k); /*!<@brief classical Gram-Schmidt orthogonalization*/
    ta[k] += V[k]*f;      /*!<@brief DGKS correlction*/
    tb[k-1] += V[k-1]*f;

    if(k==fin-1) break;
    V[k+1] = f;
  }
  f *= tb[fin-1];
}

// classical Gram-Schmidt orthogonalization
void EigenModesSolver_IRL::
orthogonalize(Field& f,const vector<Field>& V,int N)const{
  size_t size = f.size();
  assert(size%2 ==0);

  std::slice re(0,size/2,2);
  std::slice im(1,size/2,2);

  vector<double> sr(N);
  vector<double> si(N);

  for(int j=0; j<N; ++j){
    sr[j] = V[j]*f;
    si[j] = V[j].im_prod(f);
  }
  for(int j=0; j<N; ++j){
    f.add(re, -sr[j]*V[j][re] +si[j]*V[j][im]);
    f.add(im, -sr[j]*V[j][im] -si[j]*V[j][re]);
  }
}

void EigenModesSolver_IRL::setUnit(vector<double>& Qt)const{
  for(int i=0; i<Nm_*Nm_; ++i) Qt[i] = 0.0;
  for(int k=0; k<Nm_; ++k) Qt[k*Nm_+k] = 1.0;
}

int EigenModesSolver_IRL::diagonalize(vector<double>& ta,vector<double>& tb,
				       vector<double>& Qt,int Nrange)const{
  setUnit(Qt);
  int Niter = 100*Nm_;
  int kmin = 0, kmax = Nrange-1;

  for(int iter=0; iter<Niter; ++iter){
    double sft = 0.5*(ta[kmax]+ta[kmax-1]);
    double dif = 0.5*(ta[kmax]-ta[kmax-1]);
    sft += dif/fabs(dif)*sqrt(dif*dif +tb[kmax-1]*tb[kmax-1]);

    QRfact_Givens(ta,tb,Qt,Nrange,sft,kmin,kmax); // transformation

    for(int j=kmax; j>kmin; --j){   
      double dds = fabs(ta[j-1])+fabs(ta[j]);
      if(fabs(tb[j-1])+dds > dds){
        kmax = j;
	break;
      }
      if(j==kmin+1) return iter;  /*!< Convergence criterion */
    }
  }
  CCIO::cout<<"[QR method] reached to maximum iterations: "<<Niter<<"\n";
  CCIO::cout<<"failed at kmax="<<kmax<<"\n";
  abort();
}

void EigenModesSolver_IRL::QRfact_Givens(vector<double>& ta,vector<double>& tb,
					 vector<double>& Qt,int Nrange,
					 double sft,int k_min,int k_max)const{
  double x;
  for(int k=k_min; k<k_max; ++k){ // k_max is the maximum idx
    
    double r = ((k==k_min)? (ta[k]-sft)/tb[k] : tb[k-1]/x);        
    double c = r/sqrt(1.0+r*r);   //  cos_x
    double s = 1.0/sqrt(1.0+r*r); //  sin_x

    for(int i=0; i<Nrange; ++i){        
      double Qt_k  = Qt[    k*Nm_+i];
      double Qt_kp = Qt[(k+1)*Nm_+i];
      Qt[    k*Nm_+i] =  c*Qt_k +s*Qt_kp;
      Qt[(k+1)*Nm_+i] = -s*Qt_k +c*Qt_kp;
    }
    double a_k = ta[k]; 
    double b_k = tb[k]; 

    // unitary transformation G(k+1,k)^T A G(k+1,k)
    ta[k  ]= c*c*a_k +s*s*ta[k+1] +2.0*c*s*b_k; // k,k
    tb[k  ]= (c*c-s*s)*b_k +c*s*(ta[k+1] -a_k); // k,k+1
    ta[k+1]= s*s*a_k +c*c*ta[k+1] -2.0*c*s*b_k; // k+1,k+1
    
    if(k != k_min) tb[k-1]= c*tb[k-1] +s*x;     // k-1,k   
    if(k == k_max-1) break;
    x = s*tb[k+1];
    tb[k+1]*= c;                                // k+1,k+2
  }
}
