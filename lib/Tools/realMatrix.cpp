/*!@file realMatrix.cpp
  @brief implementation of the utility functions for real square matrix
*/
#include "realMatrix.hpp"
#include <vector>
#include "Communicator/comm_io.hpp"
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_permutation.h>
#include<gsl/gsl_linalg.h>

using namespace std;

namespace RealMatrix{
  void transpose(std::valarray<double>& Mt,const std::valarray<double>& M){
    int Dim = sqrt(double(M.size()));
    
    for(int i=0; i<Dim; ++i){
      Mt[Dim*i+i] = M[Dim*i+i];

      for(int j=i+1; j<Dim; ++j){
	int ij = Dim*i+j;
	int ji = Dim*j+i;

	double tmp = M[ij];
	Mt[ij] = M[ji];
	Mt[ji] = tmp;
      }
    }
  }

  void trace(double& tr,const valarray<double>& M){
    int Dim = sqrt(double(M.size()));
    tr = 0.0;
    for(int i=0; i<Dim; ++i) tr+= M[Dim*i+i];    
  }
  /*
  void invert(valarray<double>& Mi,const valarray<double>& M){
    int Dim = sqrt(double(M.size()));
    
    vector<int> ip(Dim);
    for(int i=0; i<Dim; ++i) ip[i] = i;

    valarray<double> Mt = M;

    valarray<double> vv(Dim);

    for(int i=0; i<Dim; ++i){
      double big = 0.0;
      for(int j=0; j<Dim; ++j) big = max(big,fabs(Mt[i*Dim+j]));
      vv[i] = 1.0/big;
    }

    for(int k=0; k<Dim; ++k){ 
      /// setting pivot     
      int l = k;    
      double al = vv[ip[l]]*fabs(Mt[ip[l]*Dim+k]);
      for(int i=k+1; i<Dim; ++i){
	double tmp = vv[ip[i]]*fabs(Mt[ip[i]*Dim+k]);
	if(tmp>al){ al= tmp; l= i; }
      }
      if(l != k) swap(ip[k],ip[l]);
      
      /// Gauss elimination
      Mt[ip[k]*Dim+k] = 1.0/Mt[ip[k]*Dim+k];
      
      for(int i=k+1; i<Dim; ++i){
	Mt[ip[i]*Dim+k] *= Mt[ip[k]*Dim+k];

	for(int j=k+1; j<Dim; ++j)
	  Mt[ip[i]*Dim+j] -= Mt[ip[i]*Dim+k]*Mt[ip[k]*Dim+j];
      }
    }
    
    //// Inversion of M 
    for(int j=0; j<Dim; ++j){
      double tt;
      
      for(int i=0; i<Dim; ++i){     // forward substitution 
	tt = (ip[i] == j) ? 1.0 : 0.0;
	for(int l=0; l<i; ++l)
	  tt -= Mt[ip[i]*Dim+l]*Mi[l*Dim+j];
	Mi[i*Dim+j] = tt;
      }
      for(int i=Dim-1; i>=0; --i){ // backward substitution
	tt = Mi[i*Dim+j];
	for(int l=i+1; l<Dim; ++l)
	  tt -= Mt[ip[i]*Dim+l]*Mi[l*Dim+j];

	Mi[i*Dim+j] = tt*Mt[ip[i]*Dim+i];
      }
    }
  }
  */

  /*
  void invert(valarray<double>& Mi,const valarray<double>& M){

    int Dim = sqrt(double(M.size()));
    vector<int> indx(Dim);
    valarray<double> Mt = M;
    valarray<double> vv(Dim);
    
    //// LU decomposition
    for(int i=0; i<Dim; ++i){
      double big = 0.0;
      for(int j=0; j<Dim; ++j) big = max(big,fabs(Mt[i*Dim+j]));
      vv[i] = 1.0/big;
    }

    for(int j=0; j<Dim; ++j){

      for(int i=0; i<j; ++i){
	double sum = Mt[i*Dim+j];
	for(int k=0; k<i; ++k) sum -= Mt[i*Dim+k]*Mt[k*Dim+j];
	Mt[i*Dim+j] = sum;
      }

      double big=0.0;
      int imax;
      for(int i=j; i<Dim; ++i){
	imax = j;
	double sum = Mt[i*Dim+j];
	for(int k=0; k<j; ++k) sum -= Mt[i*Dim+k]*Mt[k*Dim+j];
	Mt[i*Dim+j] = sum;
	if(double dum= vv[i]*fabs(sum) >= big){
	  big = dum;
	  imax = i;
	}
      }
      if(j!=imax){
	for(int k=0; k<Dim; ++k) swap(Mt[imax*Dim+k],Mt[j*Dim+k]);
	vv[imax] = vv[j];
      }
      indx[j]= imax;
      if(j !=Dim-1)
	for(int i=j+1; i<Dim; ++i) Mt[i*Dim+j] /= Mt[j*Dim+j];
    }

    ///////
    for(int k=0; k<Dim; ++k){
      for(int l=0; l<Dim; ++l)
	CCIO::cout<<" LU["<<k<<","<<l<<"]="<<Mt[k*Dim+l];
      CCIO::cout<<"\n";
    }


    //// forward/backward subtractions
    for(int j=0; j<Dim; ++j){
      vector<double> col(Dim,0.0);
      col[j] = 1.0;

      int ii = 0;
      for(int i=0; i<Dim; ++i){
	double sum = col[indx[i]];
	col[indx[i]] = col[i];
	if(ii)
	  for(int k=ii; k<i; ++k) sum -= Mt[i*Dim+k]*col[k];
	else if(sum) ii = i;
	col[i] = sum;
      }
      for(int i=Dim-1; i>=0; --i){
	double sum = col[i];
	for(int k=i+1; k<Dim; ++k) sum -= Mt[i*Dim+k]*col[k];
	col[i] = sum/Mt[i*Dim+i];
      }
      for(int i=0; i<Dim; ++i) Mi[i*Dim+j] = col[i];
    }
  }
  */
  void invert(valarray<double>& Mi,const valarray<double>& M){
    int Dim = sqrt(double(M.size()));
    valarray<double> Mt(M);
    gsl_matrix_view m = gsl_matrix_view_array(&(Mt[0]),Dim,Dim);
    gsl_matrix_view mi = gsl_matrix_view_array(&(Mi[0]),Dim,Dim);
    gsl_permutation* p = gsl_permutation_alloc(Dim);
    int sgn;
    
    gsl_linalg_LU_decomp(&m.matrix,p,&sgn);
    gsl_linalg_LU_invert(&m.matrix,p,&mi.matrix);
  }    
}



