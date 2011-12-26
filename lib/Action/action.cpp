/*!
  @file action.cpp
  
  @brief Definition of non virtual functions in abstract Action class

 */

#include "action.h"
#include "Communicator/comm_io.hpp"
#include "include/field.h"

void Action::monitor_force(Field& force, std::string ActionName){
  double f_re = force.average_real();
  double f_im = force.average_imag();
  double f_max= force.max_element();

  CCIO::cout<<"["<<ActionName<<"]\n";  
  //we want to check the real part of iP only for debugging purposes
  //so I'm putting the output to a higher level of verbosity (=4)
  CCIO::cout<<" +------- average MD-force (iPdot) = "  
#if VERBOSITY>3  
	    <<"(" << f_re <<","
#endif
	    << f_im 
#if VERBOSITY>3  
	    <<")"
#endif
	    << "\n";
  
  CCIO::cout<<" +------- maximum MD-force (iPdot) = "
	    << f_max
	    << "\n";

}
