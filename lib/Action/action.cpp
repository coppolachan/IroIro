/*!
  @file action.cpp
  
  @brief Definition of non virtual functions in abstract Action class

  Time-stamp: <2013-04-16 11:20:31 neo>
 */

#include "action.hpp"
#include "Communicator/comm_io.hpp"


void Action::monitor_force(GaugeField& force, std::string ActionName){
  // calculates the average of the force field
  double f_abs = force.data.average_abs(); 
  // calculates the maximum of the force field
  double f_max = force.data.max_element();

  // Print out the result 
  CCIO::cout<<"    ["<<ActionName<<"]\n";  
  CCIO::cout<<"    :::::::::: Average |MD-force| = "  
	    << f_abs
	    << "\n";
  
  CCIO::cout<<"    :::::::::: Maximum  MD-force  = "
	    << f_max
	    << "\n";
}
