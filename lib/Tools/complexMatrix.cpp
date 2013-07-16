/*!@file complexMatrix.cpp
  @brief implementation of the utility functions for complex square matrix
*/
#include "complexMatrix.hpp"
#include <complex.h>
#include <vector>

using namespace std;

namespace ComplexMatrix{
  void hermite(std::valarray<double>& Md,const std::valarray<double>& M){
    int Dim = sqrt(M.size()/2);
    
    for(int i=0; i<Dim; ++i){
      Md[2*(Dim*i+i)+1] = -M[2*(Dim*i+i)+1];

      for(int j=i+1; j<Dim; ++j){
	int ij = 2*(Dim*i+j);
	int ji = 2*(Dim*j+i);

	double re = M[ij];
	double im = M[ij+1];
	
	Md[ij  ] = M[ji  ];
	Md[ij+1] =-M[ji+1];
	
	Md[ji  ] = re;
	Md[ji+1] =-im;
      }
    }
  }

  void transpose(std::valarray<double>& Mt,const std::valarray<double>& M){
    int Dim = sqrt(M.size()/2);
    
    for(int i=0; i<Dim; ++i){
      Mt[2*(Dim*i+i)+1] = M[2*(Dim*i+i)+1];

      for(int j=i+1; j<Dim; ++j){
	int ij = 2*(Dim*i+j);
	int ji = 2*(Dim*j+i);

	double re = M[ij];
	double im = M[ij+1];

	Mt[ij  ] = M[ji  ];
	Mt[ij+1] = M[ji+1];
	
	Mt[ji  ] = re;
	Mt[ji+1] = im;
      }
    }
  }

  void conjugate(std::valarray<double>& Mc,const std::valarray<double>& M){
    int Dim = sqrt(M.size()/2);
    
    for(int i=0; i<Dim; ++i){
      Mc[2*(Dim*i+i)+1] = -M[2*(Dim*i+i)+1];

      for(int j=i+1; j<Dim; ++j){
	int ij = 2*(Dim*i+j);
	int ji = 2*(Dim*j+i);

	Mc[ij] = M[ij];
	Mc[ji] = M[ji];

	Mc[ij+1] = -M[ij+1];
	Mc[ji+1] = -M[ji+1];
      }
    }
  }

  void trace(double& tr,double& ti,const std::valarray<double>& M){
    int Dim = sqrt(M.size()/2);
    
    tr = 0.0;
    ti = 0.0;

    for(int i=0; i<Dim; ++i){
      tr += M[2*(Dim*i+i)];    
      ti += M[2*(Dim*i+i)+1];    
    }
  }

  void invert(std::valarray<double>& Mi,const std::valarray<double>& M){
    int Dim = sqrt(M.size()/2);
    
    vector<int> ip(Dim);
    vector<double> vv(Dim);
  
    valarray<double> Mt(M);
    double _Complex* Mtp = reinterpret_cast<double _Complex*>(&Mt[0]);
    double _Complex* Mip = reinterpret_cast<double _Complex*>(&Mi[0]);

    for(int i=0; i<Dim; ++i) ip[i] = i;

    for(int i=0; i<Dim; ++i){
      double big = 0.0;
      for(int j=0; j<Dim; ++j)
	big = max(big,cabs(Mtp[i*Dim+j]));
      vv[i] = 1.0/big;
    }
    
    for(int k=0; k<Dim; ++k){ 
      /// setting pivot     
      int l = k;    
      double al = vv[ip[l]]*cabs(Mtp[ip[l]*Dim+k]);
      
      for(int i=k+1; i<Dim; ++i){
	double tmp = vv[ip[i]]*cabs(Mtp[ip[i]*Dim+k]);
	if(tmp>al){ al= tmp; l= i; }
      }
      if(l != k) swap(ip[k],ip[l]);
      
      /// Gauss elimination
      Mtp[ip[k]*Dim+k] = 1.0/Mtp[ip[k]*Dim+k];
      
      for(int i=k+1; i<Dim; ++i){
	Mtp[ip[i]*Dim+k] *= Mtp[ip[k]*Dim+k];

	for(int j=k+1; j<Dim; ++j)
	  Mtp[ip[i]*Dim+j] -= Mtp[ip[i]*Dim+k]*Mtp[ip[k]*Dim+j];
      }
    }
    
    /******** Inversion of M ************/
    for(int j=0; j<Dim; ++j){
      double _Complex tt;
      
      for(int i=0; i<Dim; ++i){     // forward substitution 
	tt = (ip[i] == j) ? 1.0 : 0.0;
	for(int l=0; l<i; ++l)
	  tt -= Mtp[ip[i]*Dim+l]*Mip[l*Dim+j];
	Mip[i*Dim+j] = tt; 
      }
      for(int i=Dim-1; i>=0; --i){ // backward substitution
	tt = Mip[i*Dim+j];
	for(int l=i+1; l<Dim; ++l)
	  tt -= Mtp[ip[i]*Dim+l]*Mip[l*Dim+j];
	
	Mip[i*Dim+j] = tt*Mtp[ip[i]*Dim+i];
      }
    }
  }
}



