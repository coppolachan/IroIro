/*!
 * @file solver.cpp
 *
 * @brief Definition of SolverOutput structure
 *
 * Time-stamp: <2013-04-23 16:17:22 neo>
 */

#include "solver.hpp"
#include "Communicator/comm_io.hpp"

void SolverOutput::print(std::string FurtherMsg) {
  CCIO::cout << "Solver ["<<Msg<<"]: Iterations = "<<Iterations
	     << " Residual = "<<diff
	     << " Timing(ms) = "<< timing 
	     << "\n"<< FurtherMsg << std::endl;
}
