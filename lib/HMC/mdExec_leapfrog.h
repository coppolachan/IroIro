/*
 * @file mdExec_leapfrog.h
 *
 * @brief Declarations of MDexec_leapfrog class and Parameters 
 */
#ifndef MD_LEAPFROG_INCLUDED
#define MD_LEAPFROG_INCLUDED

#include<vector>

#include "mdExec.h"
#include "Action/action.h"
#include "Measurements/GaugeM/staples.h"
#include "include/pugi_interface.h"

struct MDexec_leapfrogParams{
  int Nexp;
  int MDsteps;
  double step_size;

  MDexec_leapfrogParams(XML::node node){
    XML::read(node, "MDsteps", MDsteps);
    XML::read(node, "step_size" , step_size);
    XML::read(node, "exp_approx", Nexp);
  }

  MDexec_leapfrogParams(int Nexp_,int MDsteps_,double step)
    :Nexp(Nexp_),MDsteps(MDsteps_),step_size(step){}
};

class MDexec_leapfrog : public MDexec {
private:
  const MDexec_leapfrogParams Params;
  const std::vector<int> Nrel_;
  const ActionSet as_;
  const Format::Format_G& gf_;
  GaugeObservers ObserverSet;
  Field* const U_;
  Field P_;

  void update_P(int lv,double ep);
  void update_U(double ep);
  void integrator_step(int level,std::vector<int>& clock);

  // Observers controls
  void register_observers();
  void attach_observer(Observer*);
  void detach_observer(Observer*){};
  void notify_observers(); 

public:
  MDexec_leapfrog(int Nexp, int MDiter, double step,
		  const ActionSet as,
		  const std::vector<int> multipliers,
		  const Format::Format_G& gf,
		  Field* const CommonF)
    :as_(as),gf_(gf),
     Params(MDexec_leapfrogParams(Nexp,MDiter,step)),
     Nrel_(multipliers),
     U_(CommonF),P_(CommonF->size())
  {
    register_observers();
  }
  
  MDexec_leapfrog(XML::node node,
		  const ActionSet as,
		  const std::vector<int> multipliers,
		  const Format::Format_G& gf,
		  Field* const CommonF)
    :as_(as),gf_(gf),
     Params(MDexec_leapfrogParams(node)),
     Nrel_(multipliers),
     U_(CommonF),P_(CommonF->size())
  {
    register_observers();
  }
  
  void init(std::vector<int>& clock,const Field& U,const RandNum& rand);
  void integrator(int level,std::vector<int>& clock);
  double calc_H() const;
  const Field get_U() const;
};

#endif  //MD_LEAPFROG_INCLUDED
