/*!
 * @file action_gaugetype_factory.cpp 
 *
 * @brief Definition of methods for Gauge-type action factories
 *
 * Time-stamp: <2013-04-17 15:15:43 neo>
 */

#include "action_gaugetype_factory.hpp"

//////////////////////////////////////////////////////////////////////////////////
WilsonGaugeActionFactory::WilsonGaugeActionFactory(const XML::node node):
  Action_node(node){};

ActionGaugeWilson* WilsonGaugeActionFactory::getGaugeAction(GaugeField* const G,
							    SmartConf* const SC)
{
  return new ActionGaugeWilson(Action_node,G); 
}
//////////////////////////////////////////////////////////////////////////////////
WilsonGaugeAdjointActionFactory::WilsonGaugeAdjointActionFactory(const XML::node node):
  Action_node(node){};

ActionGaugeWilsonAdjoint* WilsonGaugeAdjointActionFactory::getGaugeAction(GaugeField* const G, 
									  SmartConf* const SC)
{
  return new ActionGaugeWilsonAdjoint(Action_node,G); 
}
//////////////////////////////////////////////////////////////////////////////////
RectGaugeActionFactory::RectGaugeActionFactory(const XML::node node):
  Action_node(node){};

ActionGaugeRect* RectGaugeActionFactory::getGaugeAction(GaugeField* const G, 
							SmartConf* const SC)
{
  return new ActionGaugeRect(Action_node,G); 
}
//////////////////////////////////////////////////////////////////////////////////
IwasakiGaugeActionFactory::IwasakiGaugeActionFactory(const XML::node node):
  Action_node(node){};

ActionGaugeRect* IwasakiGaugeActionFactory::getGaugeAction(GaugeField* const G,
							   SmartConf* const SC){
  return new ActionGaugeRect(Action_node,
			     3.648,
			     -0.331,
			     G);
}
//////////////////////////////////////////////////////////////////////////////////
SymanzikGaugeActionFactory::SymanzikGaugeActionFactory(const XML::node node):
  Action_node(node){};

ActionGaugeRect* SymanzikGaugeActionFactory::getGaugeAction(GaugeField* const G,
							    SmartConf* const SC){
  return new ActionGaugeRect(Action_node,
			     5.0/3.0,
			     -1.0/12.0,
			     G); 
}
//////////////////////////////////////////////////////////////////////////////////
DBW2GaugeActionFactory::DBW2GaugeActionFactory(const XML::node node):
Action_node(node){};

ActionGaugeRect* DBW2GaugeActionFactory::getGaugeAction(GaugeField* const G,
							SmartConf* const SC){
  return new ActionGaugeRect(Action_node,
			     12.2704,
			     -1.4088,
			     G);
}
//////////////////////////////////////////////////////////////////////////////////
