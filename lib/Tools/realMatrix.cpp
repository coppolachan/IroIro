/*!@file realMatrix.cpp
  @brief implementation of the utility functions for real square matrix
*/
#include "realMatrix.hpp"
#include <vector>

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

  void invert(valarray<double>& Mi,const valarray<double>& M){
    int Dim = sqrt(double(M.size()));
    
    vector<int> ip(Dim);
    for(int i=0; i<Dim; ++i) ip[i] = i;

    vector<double> Mt = M;
    vector<double> vv(Dim);

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
    
    /******** Inversion of M ************/
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
}



