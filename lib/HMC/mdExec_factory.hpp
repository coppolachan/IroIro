/*!
  @file mdExec_factory.hpp
  @brief Defines the Factories for Molecular Dynamics Integrator
*/

#ifndef MDEXEC_FACT_HPP_
#define MDEXEC_FACT_HPP_

#include <iostream>
#include <string>
#include <map>
#include <vector>

#include "include/pugi_interface.h"
#include "HMC/mdExec_factory_abs.hpp"
#include "Communicator/comm_io.hpp"
#include "HMC/mdExec_leapfrog.hpp"
#include "HMC/mdExec_2MN.hpp"
#include "Smearing/smearingFactories.hpp"
/*
 * :::::::::::::::::::::   Concrete classes
 */
class MDIntegrator_LeapfrogFactory : public MDIntegratorFactory {
  ActionSetFactory ActSetFactory;
  const XML::node Integrator_node;
  SmartConfFactory SCFactory;
  GaugeField* CommonField;

  MDexec* createLeapfrog(){
    try{
      return new MDexec_leapfrog(Integrator_node,
				 ActSetFactory.
				 getActionSet(CommonField,
					      SCFactory.getSmartConfiguration()), 
				 ActSetFactory.getMultipliers(),
				 SCFactory.getSmartConfiguration());
    }catch(...){
      std::cerr << "Error in creating leapfrog integrator" << std::endl;
      abort();
    }
  }
public:
  MDexec* getMDIntegrator(){return createLeapfrog();}
  
  ~MDIntegrator_LeapfrogFactory(){delete CommonField;}

  MDIntegrator_LeapfrogFactory(XML::node node)
    :ActSetFactory(node),
     Integrator_node(node),
     SCFactory(node),
     CommonField(SCFactory.getSmartConfiguration()->ThinLinks){  } 
};

class MDIntegrator_2MN_Factory : public MDIntegratorFactory {
  ActionSetFactory ActSetFactory;
  const XML::node Integrator_node;
  SmartConfFactory SCFactory;
  GaugeField* CommonField;

  MDexec* create_2MN(){
    try{
      return new MDexec_2MN(Integrator_node,
			    ActSetFactory.
			    getActionSet(CommonField,
					 SCFactory.getSmartConfiguration()), 
			    ActSetFactory.getMultipliers(),
			    SCFactory.getSmartConfiguration());
    }catch(...){
      std::cerr << "Error in creating 2MN integrator" << std::endl;
      abort();
    }
  }
public:
  MDexec* getMDIntegrator(){return create_2MN();}
  
  ~MDIntegrator_2MN_Factory(){delete CommonField;}

  MDIntegrator_2MN_Factory(XML::node node)
    :ActSetFactory(node),
     Integrator_node(node),
     SCFactory(node),
     CommonField(SCFactory.getSmartConfiguration()->ThinLinks){  } 
};

#endif // MDEXEC_FACT_HPP_
