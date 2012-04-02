/*! 
 * @file action_Nf2_DomainWall.hpp
 * @brief Declaration of Action_Nf2_DomainWall class
 */
#ifndef ACTION_NF2_DOMAINWALL_INCLUDED
#define ACTION_NF2_DOMAINWALL_INCLUDED

#include "Action/action_Nf2_ratio.hpp"
#include "Dirac_ops/dirac_DomainWall.hpp"

class Action_Nf2_DomainWall :public Action{
private:
  Action_Nf2_ratio action_;
public:
  Action_Nf2_DomainWall(GaugeField* const F,
			Dirac_optimalDomainWall* const D,
			Dirac_optimalDomainWall* const Dpv,
			const Solver* Solv,
			const Solver* SolvPv,
			SmartConf* smart_conf = NULL)
    :action_(F,D,Dpv,Solv,SolvPv,smart_conf){}
  
  void observer_update();
  void init(const RandNum& rand);

  double calc_H();
  GaugeField md_force();
};
#endif
