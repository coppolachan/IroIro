/*!
  @file mdExec_Factory.h
  
  @brief Defines the Factories for Molecular Dynamics Integrator
 
*/

#ifndef MDEXEC_FACT_H_
#define MDEXEC_FACT_H_

#include <iostream>
#include <string>
#include <map>
#include <vector>

#include "include/pugi_interface.h"
#include "Tools/RAIIFactory.hpp"

#include "Action/action_Factory.h"
#include "Action/action_gaugetype_factory.hpp"
#include "Action/action_fermiontype_factory.hpp"
#include "Communicator/comm_io.hpp"

#include "HMC/mdExec_leapfrog.h"

/*!
 *@class IntegratorFactory
 *
 *@brief Abstract Factory class for Molecular Dynamics Integrator
 *
 */
class MDIntegratorFactory {
public:
  /*!
   Virtual function returning the Molecular Dynamics Integrator
   */
  virtual MDexec* getMDIntegrator() = 0;
};

/*!
 *@class ActionSetFactory
 *
 *@brief Factory class for grouping Action in levels
 *
 */
class ActionSetFactory {

  std::multimap< int, ActionFactory* > LevelMap;/*<- Map to store the factory
						 * levels
						 */
  std::vector<int> multipliers; 

  RaiiFactoryObjVec<Action> ActionPointers;

  /*!
    Function returning the ActionSet
  */
  void next_step(XML::node node) {
    int m;
    XML::descend(node, "step");
    if (node !=NULL) {
      
      XML::read(node,"multiplier", m);
      multipliers.push_back(m);
      
      CCIO::cout << "Level : " << multipliers.size() 
		 << " with multiplier " << m << std::endl;
      
      for (XML::node action = node.child("Action"); 
	   action; 
	   action = action.next_sibling("Action")) 
	{
	  std::string type =action.attribute("type").value() ;
	  std::string action_name =action.attribute("name").value() ;
	  CCIO::cout << "Found Action  [" << action_name
		     << "] of type [" << type << "]" << std::endl;
	  
	  //Associate name to the correct action creator
	  if(!type.compare("Gauge")){
	    LevelMap.insert( std::pair<int,ActionFactory*>
			     (multipliers.size(), 
			      GaugeAction::createGaugeActionFactory(action) ));
	      }

	    if(!type.compare("Fermion")){
	      LevelMap.insert( std::pair<int,ActionFactory*>
			     (multipliers.size(), 
			      FermionAction::createFermionActionFactory(action) ));

	    
	    }

	}
      next_step(node);
    }
    
    
  }
  
  ActionSetFactory(const ActionSetFactory&){} //hide copy constructor

public:  
  ActionSet getActionSet(const Format::Format_G& GaugeForm,
			  Field* const GaugeField){
    std::multimap<int , ActionFactory* >::iterator it;

    ActionSet Levels;
    CCIO::cout<< "Number of actions " << LevelMap.size() 
	     << "  divided in " << multipliers.size() << " nested levels " 
	     << std::endl;

    Levels.resize(multipliers.size());
 
    for ( it = LevelMap.begin() ; it != LevelMap.end(); it++) {
      Levels[(*it).first-1].
	push_back(ActionPointers.save((*it).second->getAction(GaugeForm,GaugeField)));
    }

    return Levels;
  };

  const std::vector<int> getMultipliers(){
    return multipliers;
  } 
  
  ActionSetFactory(XML::node Integrator_node) {
    //traverse the xml input looking for <step> objects
    next_step(Integrator_node);
  };

};//end of class ActionSetFactory definition



/*
 * :::::::::::::::::::::   Concrete classes
 */

class MDIntegrator_LeapfrogFactory : public MDIntegratorFactory {
  ActionSetFactory ActSetFactory;
  const Format::Format_G& Format;
  const XML::node Integrator_node;
  Field* CommonField;
  RaiiFactoryObj<Staples> StaplePointer;

  MDexec* createLeapfrog() {
    
    CommonField->resize(Format.size()); 
    try {
      StaplePointer.save(new Staples(Format));
      return new MDexec_leapfrog(Integrator_node,
				 ActSetFactory.getActionSet(Format, CommonField), 
				 ActSetFactory.getMultipliers(),
				 StaplePointer.get(),
				 Format,
				 CommonField);
      
    } catch (...) {
      std::cerr << "Error in creating leapfrog" << std::endl;
      abort();
    }
    
    
    
  };
  
public:
  MDexec* getMDIntegrator () {
    return createLeapfrog();
  };
  
  ~MDIntegrator_LeapfrogFactory() {
    delete CommonField;
  }

  MDIntegrator_LeapfrogFactory (XML::node node, Format::Format_G& F):
    ActSetFactory(node),
    Integrator_node(node),
    Format(F),
    CommonField(new Field){
  } 

};

namespace Integrators{
  static MDIntegratorFactory* Integr;
  MDIntegratorFactory* createIntegratorFactory(XML::node,
					       Format::Format_G&);
  
}

#endif // MDEXEC_FACT_H_
