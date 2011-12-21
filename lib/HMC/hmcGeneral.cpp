//--------------------------------------------------------------------
// hmcGeneral.cpp
//--------------------------------------------------------------------
#include "hmcGeneral.h"
#include "mdExec.h"

#include "include/field.h"
#include "include/format_G.h"

#include "Tools/sunMat.h"
#include "Communicator/comm_io.hpp"
#include "Action/action.h"

using namespace std;

void HMCgeneral::evolve(Field& Uin)const{

  for(int iter=1; iter <= Params.ThermalizationSteps; ++iter){
    CCIO::cout << "-- # Thermalization step = "<< iter << endl;

    double Hdiff = evolve_step(Uin);
    CCIO::cout<< "dH = "<<Hdiff << std::endl;
    Uin = md_->get_U();  //accept every time
  }

  for(int iter=1; iter <= Params.Nsweeps; ++iter){
    CCIO::cout << "-- # Sweep = "<< iter << endl;

    double Hdiff = evolve_step(Uin);
    
    if(metropolis_test(Hdiff)) Uin = md_->get_U();
  }
}

double HMCgeneral::evolve_step(Field& Uin)const{
    
    vector<int> clock;
    md_->init(clock,Uin,*rand_);     // set U and initialize P and phi's 

    int init_level = 0;
    double H0 = md_->calc_H();     // current state            
    CCIO::cout<<"total H_before = "<< H0 <<endl;

    md_->integrator(init_level,clock);

    double H1 = md_->calc_H();     // updated state            
    CCIO::cout<<"Total H_after = "<< H1 <<endl;

    return (H1-H0);
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

