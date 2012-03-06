/*! 
 * @file action_Nf2_DomainWall.hpp
 *
 * @brief Declaration of Action_Nf2_DomainWall class
 *
 */
#ifndef ACTION_NF2_DOMAINWALL_INCLUDED
#define ACTION_NF2_DOMAINWALL_INCLUDED

//#include "Tools/randNum_MP.h"
#include "Action/action_Nf2_ratio.hpp"
#include "Dirac_ops/dirac_DomainWall.hpp"

/*!
 * @class Action_Nf2_DomainWall
 *
 * @brief Class to calculate TwoFlavors HMC action term for DWF
 *
 */
class Action_Nf2_DomainWall :public Action{
private:
  Action_Nf2_ratio action_;
public:
  /*!
   * @brief Standard constructor 
   */
  Action_Nf2_DomainWall(GaugeField* const F,
			const Dirac_optimalDomainWall* D,
			const Dirac_optimalDomainWall* Dpv,
			const Solver* Solv,
			const Solver* SolvPv)
    :action_(F,D,Dpv,Solv,SolvPv){}
  
  void observer_update(){};
  void init(const RandNum& rand);

  double calc_H();
  GaugeField md_force();
};
#endif
