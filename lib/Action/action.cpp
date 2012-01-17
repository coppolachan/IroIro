/*!
  @file action.cpp
  
  @brief Definition of non virtual functions in abstract Action class

 */

#include "action.hpp"
#include "Communicator/comm_io.hpp"
#include "include/field.h"

void Action::monitor_force(Field& force, std::string ActionName){
  double f_abs = force.average_abs();
  double f_max = force.max_element();

  CCIO::cout<<"    ["<<ActionName<<"]\n";  
  //we want to check the real part of iP only for debugging purposes
  CCIO::cout<<" +------- average |MD-force| (iPdot) = "  
	    << f_abs
	    << "\n";
  
  CCIO::cout<<" +------- maximum  MD-force  (iPdot) = "
	    << f_max
	    << "\n";
}
