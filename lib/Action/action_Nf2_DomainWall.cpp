/*! 
  @file action_Nf2_DomainWall.cpp
  @brief Definition of Action_Nf2_DomainWall class
  
  Time-stamp: <2013-04-22 16:59:42 neo>
 */
#include "action_Nf2_DomainWall.hpp"

void Action_Nf2_DomainWall::init(const RandNum& rand){
  action_.init(rand);
}
double Action_Nf2_DomainWall::calc_H(){
  return action_.calc_H();
}
void Action_Nf2_DomainWall::observer_update(){
  action_.observer_update();
}

GaugeField Action_Nf2_DomainWall::md_force(){
  return action_.md_force();
}
