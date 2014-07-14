#ifndef SIGNAPPROX_CHEBYSHEV_INCLUDED
#define SIGNAPPROX_CHEBYSHEV_INCLUDED

#include <vector>
#include <valarray>

class SignApprox_Chebyshev{
  int Npoly_,N_;
  double eps_;

  double precYs_,stpCnd_;
  int Nsmpl_,maxIter_,maxNewton_;

  double u_;    
  std::vector<double> c_; /// coefficients of the Chebyshef series
  std::vector<double> y_; /// extrema of the polynomial

  void doApprox();
  void updateCs(int it);
  void initYs();
  void updateYs(int it);
  void setCbMat(std::valarray<double>&);

  double newtonStep(int,double,
		    const std::vector<double>&,
		    const std::vector<double>&)const;
  double maxh()const;

  double y2z(double y)const{ return (2.0*y -1.0-eps_)/(1.0-eps_);}
  double z2y(double z)const{ return 0.5*((1.0 -eps_)*z +1.0 +eps_);}

public:
  SignApprox_Chebyshev(int Npoly,double epsilon,
		       double precYs=1.0e-5,
		       double stpCnd=1.0e-2,
		       int Nsmpl = 100,
		       int maxIter=40,
		       int maxNewton = 50)
    :Npoly_(Npoly),N_(Npoly+2),eps_(epsilon),
     precYs_(precYs),stpCnd_(stpCnd),Nsmpl_(Nsmpl),
     maxIter_(maxIter),maxNewton_(maxNewton),
     y_(N_),c_(Npoly+1){  doApprox(); }

  std::vector<double> getCoeffs()const{ return c_;}
  double func(double y)const;
};

#endif
