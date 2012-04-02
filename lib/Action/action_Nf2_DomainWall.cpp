/*! 
 * @file action_Nf2_DomainWall.cpp
 * @brief Definition of Action_Nf2_DomainWall class
 */
#include "action_Nf2_DomainWall.hpp"
#include "include/messages_macros.hpp"

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
