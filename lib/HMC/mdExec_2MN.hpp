/*
 * @file mdExec_2MN.hpp
 * @brief Declarations of MDexec_2MN class and Parameters
 * 
 * It performs the 2nd order minimum norm integration  
 */
#ifndef MD_2MN_INCLUDED
#define MD_2MN_INCLUDED

#include "mdExec.hpp"
#include "Action/action.hpp"
#include "include/pugi_interface.h"
#include "Smearing/smartConf.hpp"
#include<vector>

struct MDexec_2MN_Params{
  int Nexp;
  int MDsteps;
  double step_size;
  double lambda;



  MDexec_2MN_Params(XML::node node){
    XML::read(node, "MDsteps", MDsteps);
    XML::read(node, "step_size" , step_size);
    XML::read(node, "exp_approx", Nexp);
    lambda = 0.1931833275037836;
  }

  MDexec_2MN_Params(int Nexp_,int MDsteps_,double step)
    :Nexp(Nexp_),MDsteps(MDsteps_),step_size(step),
     lambda(0.1931833275037836){}
};

class MDexec_2MN : public MDexec {
private:
  const MDexec_2MN_Params Params;
  const std::vector<int> Nrel_;
  const ActionSet as_;
  ObserverList observers_;
  GaugeField* const U_;
  GaugeField P_;

  // Private functions
  void update_P(int lv,double ep);
  void update_U(double ep);
  void integrator_step(int level,std::vector<int>& clock);

  // Observers controls
  void register_observers();
  void notify_observers(); 

public:
  MDexec_2MN(int Nexp, int MDiter, double step,
	     const ActionSet as,
	     const std::vector<int> multipliers,
	     SmartConf* const CommonU)
    :as_(as),
     Params(MDexec_2MN_Params(Nexp,MDiter,step)),
     Nrel_(multipliers),
     U_(CommonU->ThinLinks){
    observers_.push_back(CommonU);//Attach smearing as 1st observer
    register_observers();
  }
  
  MDexec_2MN(int Nexp, int MDiter, double step,
	     const ActionSet as,
	     const std::vector<int> multipliers,
	     GaugeField* const CommonU)
    :as_(as),
     Params(MDexec_2MN_Params(Nexp,MDiter,step)),
     Nrel_(multipliers),
     U_(CommonU),P_(){
    register_observers();
  }
  
  MDexec_2MN(XML::node node,
	     const ActionSet as,
	     const std::vector<int> multipliers,
	     GaugeField* const CommonU)
    :as_(as),
     Params(MDexec_2MN_Params(node)),
     Nrel_(multipliers),
     U_(CommonU),P_(){
    register_observers();
  }
  
  MDexec_2MN(XML::node node,
	     const ActionSet as,
	     const std::vector<int> multipliers,
	     SmartConf* const CommonU)
    :as_(as),
     Params(MDexec_2MN_Params(node)),
     Nrel_(multipliers),
     U_(CommonU->ThinLinks){
    observers_.push_back(CommonU);//Attach smearing as 1st observer
    register_observers();
  }
  
  
  void init(std::vector<int>& clock,const GaugeField& U,const RandNum& rand);
  void integrator(int level,std::vector<int>& clock);
  double calc_H() const;
  const GaugeField get_U() const;
};

#endif  //MD_2MN_INCLUDED
