#ifndef SIGNAPPROX_CHEBYSHEV_INCLUDED
#define SIGNAPPROX_CHEBYSHEV_INCLUDED

#include <vector>

class SignApprox_Chebyshev{
  int Npoly_,N_;
  double eps_;

  double precYs_,stpCnd_;
  int Nsmpl_,maxIter_;

  double u_;    
  std::valarray<double> c_; /// coefficients of the Chebyshef series
  std::valarray<double> y_; /// extrema of the polynomial

  void doApprox();
  void updateCs(int it,const vector<double>& y);
  void updateYs(std::vector<double>& y1,const vector<double>& y0,int it);
  void setCbMat(std::vector<double>&,const vector<double>&);

  double newtonStep(double,
		    const vector<double>&,
		    const vector<double>&,
		    const vector<double>&)const;
  double maxh()const;
public:
  SignApprox_Chebyshev(int Npoly,double epsilon,
		       double precYs=1.0e-5,
		       double stpCnd=1.0e-2,
		       int Nsmpl = 100,
		       int maxIter=40,
		       int maxNewton = 50)
    :Npoly_(Npoly),N_(Npoly+2),eps_(epsilon),
     precYs_(precYs),stpCnd_(stpCnd),Nsmpl_(Nsmpl),maxIter_(maxIter),
     y_(0.0,N_),c_(0.0,N_-1){  doApprox(); }

  std::vector<double> getCoeffs(){ return c_;}
  double func(double y)const;
};

#endif
