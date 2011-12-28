/*!
 * @file solver.cpp
 *
 * @brief Definition of SolverOutput structure
 *
 *
 */

#include "solver.hpp"
#include "Communicator/comm_io.hpp"

void SolverOutput::print(std::string FurtherMsg) {
  CCIO::cout << "Solver ["<<Msg<<"]: Iterations = "<<Iterations
	     << " Residual = "<<diff
	     << " Timing(ms) = "<< timing << std::endl;
}
