#include "signApprox_Chebyshev.hpp"
#include "Tools/realMatrix.hpp"
#include "chebyshev.hpp"
#include "include/macros.hpp"

using namespace std;

double SignApprox_Chebyshev::func(double y)const{ 
  return 1.0-sqrt(y)*Chebyshev::series(y,c_);}

void SignApprox_Chebyshev::doApprox(){

  bool converged = false;
  for(int it=0; it< maxIter_; ++it){
    updateYs(it);
    updateCs(it);  

    mh = maxh();
    CCIO::cout<<"|u|="<<fabs(u_)<<" max|h|="<<mh<<"\n";

    if((mh-fabs(u_))/fabs(u_) < stpCnd_){
      converged = true;
      return;
    }
  }
  if(!conversed){
    CCIO::cout<<"did not converge\n";
    abort();
  }
}

double SignApprox_Chebyshev::maxh()const{
  double hm = 0.0;
  for(int i=0; i<=Nsmpl_; ++i){
    double y = eps_+double(i)*(1.0-eps_)/Nsmpl_;
    double sgn = fabs(func(y));
    if(hm < sgn) hm = sgn;
  }
  return hm;
}

void SignApprox_Chebyshev::updateCs(int it){
  vector<double> A(N_*N_), Ainv(N_*N_);
  setCbMat(A,y_);
  RealMatrix::invert(Ainv,A);
  u_= 0.0;

  CCIO::cout<<"iteration "<<it<<"\n";

  for(int i=0; i<N_; ++i) u_+= Ainv[i];
  CCIO::cout<<"u= "<<u_<<"\n";

  for(int l=0; l<N_-1; ++l){
    c_[l] = 0.0;
    for(int i=0; i<N_; ++i) c_[l]+= Ainv[N_*l+i];
  }
}

void SignApprox_Chebyshev::setCbMat(vector<double>& A){
  for(int l=0; l<y.size(); ++l){
    double z = (2*y[l]-1.0-eps_)/(1.0-eps_);
    std::vector<double> Tcb(y.size()-1); /// y.size()==Npoly_+2
    get_Chebyshev(Tcb,z);

    A[l*y.size()] = 1.0 -2.0*(l%2);
    for(int m=1; m<y.size(); ++m)
      A[l*y.size()+m] = sqrt(y[l])*Tcb[m-1];
  }
}

void SignApprox_Chebyshev::updateYs(int it){
  if(!it){
    for(int l=0; l<N_; ++l) y_[l]=cos(l*PI/(Npoly_+1));
  }else{

    vector<double> c1 = Chebyshev::dcoeff(c_);
    vector<double> c2 = Chebyshev::dcoeff(c0);

    vector<double> cd(c_.size());
    vector<double> cdd(c_.size());

    vector<double> y1;    

    for(int l=0; l<N_; ++l){
      double y0 = y_[l];

      for(int it=0; it<maxNewton_; ++it){
	
	for(int i=0; i<c_.size(); ++i) cd[i] = -0.5/sqrt(y0)*c_[i];
	for(int i=0; i<c1.size(); ++i) cd[i] -= sqrt(y0)*c1[i];

	for(int i=0; i<c_.size(); ++i) cdd[i] = 0.25/y0/sqrt(y0)*c_[i];
	for(int i=0; i<c1.size(); ++i) cdd[i] -= 1.0/sqrt(y0)*c1[i];
	for(int i=0; i<c2.size(); ++i) cdd[i] -= sqrt(y0)*c2[i];

	double yt = newtonStep(l,y0,cd,cdd);
	  
	if(fabs((yt-y0)/y0) < precYs_){
	  y_[l] = yt;
	  break;
	}else{
	  y0 = yt;
	}
      }
    }
  }
}

double SignApprox_Chebyshev::newtonStep(int Ns,double y0,
					const vector<double>& cd,
					const vector<double>& cdd)const{
  double hd = Chebyshev::series(y0,cd);
  double hdd = Chebyshev::series(y0,cdd);
  double den = 0.0;
  for(int i=0; i<Ns; ++i) den += 1.0/(y0 -y_[i]);

  return y0 -hd/(hdd-hd*den);
}

