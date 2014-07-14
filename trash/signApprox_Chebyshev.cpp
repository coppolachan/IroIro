#include "signApprox_Chebyshev.hpp"
#include "Tools/realMatrix.hpp"
#include "chebyshevUtils.hpp"
#include "Communicator/comm_io.hpp"
#include <algorithm>

using namespace std;

double SignApprox_Chebyshev::func(double y)const{ 
  return 1.0-sqrt(y)*Chebyshev::series(y,c_);}

void SignApprox_Chebyshev::doApprox(){

  bool converged = false;
  initYs();
  for(int it=0; it< maxIter_; ++it){
    updateCs(it);  
    updateYs(it);

    double mh = maxh();
    CCIO::cout<<"|u|="<<fabs(u_)<<" max|h|="<<mh<<"\n";

    if((mh-fabs(u_))/fabs(u_) < stpCnd_){
      converged = true;
      return;
    }
  }
  if(!converged){
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
  CCIO::cout<<"updateCs called\n";
  valarray<double> A(N_*N_), Ainv(N_*N_);
  setCbMat(A);
  RealMatrix::invert(Ainv,A);

  for(int k=0; k<N_; ++k){
    for(int l=0; l<N_; ++l){
      double tmp =0.0;
      for(int m=0; m<N_; ++m) tmp += A[k*N_+m]*Ainv[m*N_+l];
      CCIO::cout<<" I["<<k<<","<<l<<"]="<<tmp;
      //CCIO::cout<<" Ai["<<k<<","<<l<<"]="<<Ainv[k*N_+l];
    }
    CCIO::cout<<"\n";
  }

  u_= 0.0;

  CCIO::cout<<"iteration "<<it<<"\n";

  for(int i=0; i<N_; ++i) u_+= Ainv[i];
  CCIO::cout<<"u= "<<u_<<"\n";

  for(int l=0; l<N_-1; ++l){
    c_[l] = 0.0;
    for(int i=0; i<N_; ++i) c_[l]+= Ainv[N_*l+i];
  }
}

void SignApprox_Chebyshev::setCbMat(valarray<double>& A){
  for(int l=0; l<N_; ++l){

    std::vector<double> Tcb(c_.size()); 
    Chebyshev::func(y2z(y_[l]),Tcb);

    A[l*N_] = 1.0 -2.0*(l%2);
    //CCIO::cout<<"A["<<l<<",0]="<<A[l*N_];
    
    for(int m=0; m<c_.size(); ++m){
      A[l*N_+m+1] = sqrt(y_[l])*Tcb[m];
      //CCIO::cout<<" A["<<l<<","<<m+1<<"]="<<A[l*N_+m+1];
    }
    //CCIO::cout<<"\n";
  }
}

void SignApprox_Chebyshev::initYs(){
  for(int l=0; l<N_; ++l)
    y_[l] = z2y(Chebyshev::extrema(l,Npoly_+1));
  sort(y_.begin(),y_.end());
}

void SignApprox_Chebyshev::updateYs(int it){

  vector<double> c1 = Chebyshev::dcoeff(c_);
  vector<double> c2 = Chebyshev::dcoeff(c1);
  double dzdy= 2.0/(1.0-eps_);

  for(int i=1;i<c2.size();++i){
    double a=1.0;
    if(i==1) a = 2.0;
    CCIO::cout<<" c_["<<i<<"]="<<c_[i]
	      <<" (c1["<<i-1<<"]-c1["<<i+1<<"])/2/"<<i<<"="
	      <<  (c1[i-1]*a-c1[i+1])/2/i<<"\n";
  }
    /*
  for(int i=0;i<c2.size();++i)
  CCIO::cout<<" c_["<<i<<"]="<<c_[i]
  <<" c1["<<i<<"]="<<c1[i]
  <<" c2["<<i<<"]="<<c2[i]<<"\n";
    */
  vector<double> cd(c_.size());
  vector<double> cdd(c_.size());

  vector<double> y1;    

  for(int l=0; l<N_; ++l){
    double y0 = y_[l];

    for(int in=0; in<maxNewton_; ++in){
      
      for(int i=0; i<c_.size(); ++i) cd[i] = -0.5/sqrt(y0)*c_[i];
      for(int i=0; i<c1.size(); ++i) cd[i] -= sqrt(y0)*dzdy*c1[i];

      for(int i=0; i<c_.size(); ++i) cdd[i] = 0.25/y0/sqrt(y0)*c_[i];
      for(int i=0; i<c1.size(); ++i) cdd[i] -= 1.0/sqrt(y0)*dzdy*c1[i];
      for(int i=0; i<c2.size(); ++i) cdd[i] -= sqrt(y0)*dzdy*dzdy*c2[i];
      
      double yt = newtonStep(l,y0,cd,cdd);

      CCIO::cout<<"l="<<l<<" in="<<in<<" y0="<<y0<<"\n";

      if(fabs((yt-y0)/y0) < precYs_){
	y_[l] = yt;
	break;
      }else{
	y0 = yt;
      }
    }
  }
}

double SignApprox_Chebyshev::newtonStep(int Ns,double y0,
					const vector<double>& cd,
					const vector<double>& cdd)const{
  /*
  double h = Chebyshev::series(y2z(y0),c_);
  double hp = Chebyshev::series(y2z(y0+1.0e-6),c_);
  double hn = Chebyshev::series(y2z(y0-1.0e-6),c_);
  h *= -sqrt(y0);
  hp *= -sqrt(y0+1.0e-6);
  hn *= -sqrt(y0-1.0e-6);
  
  double hdn = (hp-h)/1.0e-6;
  double hdn2 = (hp-hn)/2.0e-6;
  */
  double hd = Chebyshev::series(y2z(y0),cd);
  double hdd = Chebyshev::series(y2z(y0),cdd);

  double den = 0.0;
  //for(int i=0; i<Ns; ++i) den += 1.0/(y0 -y_[i]);

  //  CCIO::cout<<"hdn="<<hdn<<" hdn2="<<hdn2<<" hd="<<hd<<" hdd="<<hdd<<" den="<<den<<"\n";
  return y0 -hd/(hdd-hd*den);
}

