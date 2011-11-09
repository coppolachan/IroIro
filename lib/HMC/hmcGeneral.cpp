//--------------------------------------------------------------------
// hmcGeneral.cpp
//--------------------------------------------------------------------
#include "hmcGeneral.h"
#include "mdExec.h"

#include "include/field.h"
#include "include/format_G.h"

//#include "Tools/randNum.h"
#include "Tools/sunMat.h"
#include "Communicator/comm_io.hpp"
#include "Action/action.h"


using namespace std;

void HMCgeneral::evolve(Field& Uin)const{


  for(int swp=0; swp < Params.Nsweeps; ++swp){
    CCIO::cout << "-- # sweep = "<< swp << endl;
    Field U(Uin);
    Field P(Uin.size());
    vector<int> clock;

    md_->init(clock,P,U,*rand_);    // generate internal P and phi's 

    int init_level = 0;
    double H0 = calc_H(P);            // current state            
    CCIO::cout<<"total H_before = "<< H0 <<endl;

    md_->integrator(P,init_level,clock);

    double H1 = calc_H(P);            // updated state            
    CCIO::cout<<"Total H_after = "<< H1 <<endl;
    
    if(metropolis_test(H1-H0)) Uin = U;  
  }
}

double HMCgeneral::calc_H(const Field& P)const{
  
  double pnorm = P.norm();
  CCIO::cout<<"P.norm()= "<<pnorm<<endl;
  using namespace SUNmat_utils;

  // gauge sector
  double H_local = 0.0;
  Format::Format_G gf(CommonPrms::instance()->Nvol());

  for(int site = 0; site < CommonPrms::instance()->Nvol(); ++site){
    for(int dir = 0; dir < CommonPrms::instance()->Ndim(); ++dir){
      SUNmat Pxm(P[gf.cslice(0,site, dir)]);
      H_local += ReTr(Pxm*Pxm);
    }
  }
  H_local *= -0.5;
  CCIO::cout<<"H_local="<<H_local<<endl;

  double H = Communicator::instance()->reduce_sum(H_local);
  CCIO::cout<<"H_kin = "<< H<<endl;
  
  H += md_->calc_H();
  return H;
}

bool HMCgeneral::metropolis_test(const double Hdiff)const{

  double prob = exp(-Hdiff);
  double rn = rand_->get();
  CCIO::cout<< "--------------------------------------------\n";
  CCIO::cout<< "dH = "<<Hdiff << "  Random = "<< rn 
	    << "\nAcc. Probability = " << ((prob<1.0)? prob: 1.0)<< "   ";
  if(prob >1.0){       // accepted
    CCIO::cout<<"-- ACCEPTED\n";
    return true;
  }else if(rn <= prob){// accepted
    CCIO::cout <<"-- ACCEPTED\n";
    return true;
  }else{               // rejected
    CCIO::cout <<"-- REJECTED\n";
    return false;
  }
}

