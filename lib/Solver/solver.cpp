/*!
 * @file solver.cpp
 *
 * @brief Definition of SolverOutput structure
 *
 * Time-stamp: <2013-04-22 17:40:28 neo>
 */

#include "solver.hpp"
#include "Communicator/comm_io.hpp"

void SolverOutput::print(std::string FurtherMsg) {
  CCIO::cout << "Solver ["<<Msg<<"]: Iterations = "<<Iterations
	     << " Residual = "<<diff
	     << " Timing(ms) = "<< timing << std::endl;
}
