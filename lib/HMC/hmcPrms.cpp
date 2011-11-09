//--------------------------------------------------------------------
// hmcPrms.cpp
//--------------------------------------------------------------------
#include "hmcPrms.h"
#include <iostream>
#include <cstdlib>
using namespace std;

double HMCprms::epsilon_;
int    HMCprms::MDsteps_;
int    HMCprms::Nsweeps_;
int    HMCprms::Nexp_;
double HMCprms::s_cond_;
int    HMCprms::Niter_;
std::vector<int>  HMCprms::Nrel_;

HMCprms* HMCprms::instance_;

HMCprms* HMCprms::instance(const HMC& hmc){
  if(instance_==NULL){
    instance_= new HMCprms;
    setup(hmc);
  }else{
    return instance_;
  }
}

HMCprms* HMCprms::instance(){
  if(instance_==NULL){
    cout<<"HMCprms is not initialized."<<endl;
    exit(EXIT_FAILURE);
  }
  return instance_;
}

void HMCprms::setup(const HMC& hmc){
  epsilon_= hmc.epsilon;
  MDsteps_= hmc.MDsteps;
  Nsweeps_= hmc.Nsweeps;
  Nexp_= hmc.Nexp;
  
  s_cond_= hmc.s_cond;
  Niter_= hmc.Niter;
 
  int size = hmc.Nlevel;
  for(int i = 0; i<size;++i)
    Nrel_.push_back(hmc.Nrel[i]);
}
