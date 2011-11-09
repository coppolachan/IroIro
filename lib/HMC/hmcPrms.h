//--------------------------------------------------------------------
// hmcPrms.h
//--------------------------------------------------------------------
#ifndef HMCPRMS_INCLUDED
#define HMCPRMS_INCLUDED

#include <vector>

struct HMC{
  // md step
  double epsilon;
  int   MDsteps;
  int   Nsweeps;
  int   Nexp;
  // solver
  double s_cond;
  int   Niter;
  int*  Nrel;  
  int   Nlevel;
};

class HMCprms{
private:
  // md step
  static double epsilon_;
  static int   MDsteps_;
  static int   Nsweeps_;
  static int   Nexp_;
  // solver  
  static double  s_cond_;
  static int     Niter_;
  static std::vector<int>  Nrel_;

  HMCprms(){}
  HMCprms(const HMCprms&){}
  HMCprms& operator=(const HMCprms&){}
  
  static void setup(const HMC& hmc);
  static HMCprms* instance_;

public:
  static HMCprms* instance(const HMC& hmc);
  static HMCprms* instance();

  static double epsilon(){return epsilon_;}
  static int    MDsteps(){return MDsteps_;}
  static int    Nsweeps(){return Nsweeps_;}
  static int       Nexp(){return Nexp_;}

  static double  s_cond(){return s_cond_;}
  static int      Niter(){return Niter_;}
  static std::vector<int> Nrel(){return Nrel_;}
};

#endif
