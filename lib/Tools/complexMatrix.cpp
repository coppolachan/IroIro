/*!@file complexMatrix.cpp
  @brief implementation of the utility functions for complex square matrix
*/
#include "complexMatrix.hpp"
#include <complex>
#include <vector>

using namespace std;

namespace ComplexMatrix{
  void hermite(std::valarray<double>& Md,const std::valarray<double>& M){
    int Dim = sqrt(double(M.size()/2));
    
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
    int Dim = sqrt(double(M.size()/2));
    
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
    int Dim = sqrt(double(M.size()/2));
    
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

  void trace(double& tr,double& ti,const valarray<double>& M){
    int Dim = sqrt(double(M.size()/2));
    
    tr = 0.0;
    ti = 0.0;

    for(int i=0; i<Dim; ++i){
      tr += M[2*(Dim*i+i)];    
      ti += M[2*(Dim*i+i)+1];    
    }
  }

  void invert(valarray<double>& Mi,const valarray<double>& M){
    int Dim = sqrt(double(M.size()/2));
    
    vector<int> ip(Dim);
    for(int i=0; i<Dim; ++i) ip[i] = i;

    vector<complex<double> > Mt(Dim);
    for(int i=0; i< Dim; ++i)
      Mt[i] = complex<double>(M[2*i],M[2*i+1]);
  
    vector<double> vv(Dim);
    for(int i=0; i<Dim; ++i){
      double big = 0.0;
      for(int j=0; j<Dim; ++j)
	big = max(big,abs(Mt[i*Dim+j]));
      vv[i] = 1.0/big;
    }
    
    for(int k=0; k<Dim; ++k){ 
      /// setting pivot     
      int l = k;    
      double al = vv[ip[l]]*abs(Mt[ip[l]*Dim+k]);
      
      for(int i=k+1; i<Dim; ++i){
	double tmp = vv[ip[i]]*abs(Mt[ip[i]*Dim+k]);
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
      complex<double> tt;
      
      for(int i=0; i<Dim; ++i){     // forward substitution 
	tt = (ip[i] == j) ? 1.0 : 0.0;
	for(int l=0; l<i; ++l)
	  tt -= Mt[ip[i]*Dim+l]*complex<double>(Mi[2*(l*Dim+j)],Mi[2*(l*Dim+j)+1]);

	Mi[2*(i*Dim+j)  ] = tt.real(); 
	Mi[2*(i*Dim+j)+1] = tt.imag(); 
      }
      for(int i=Dim-1; i>=0; --i){ // backward substitution
	tt = complex<double>(Mi[2*(i*Dim+j)],Mi[2*(i*Dim+j)+1]);
	for(int l=i+1; l<Dim; ++l)
	  tt -= Mt[ip[i]*Dim+l]*complex<double>(Mi[2*(l*Dim+j)],Mi[2*(l*Dim+j)+1]);

	tt *= Mt[ip[i]*Dim+i];
	
	Mi[2*(i*Dim+j)  ] = tt.real();
	Mi[2*(i*Dim+j)+1] = tt.imag();
      }
    }
  }
}



