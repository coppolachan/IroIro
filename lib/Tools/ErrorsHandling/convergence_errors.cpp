/*
 * @file convergence_errors.cpp
 *
 * @brief Collection of simple convergence error handling functions
 *
 * Time-stamp: <2013-04-23 11:10:42 neo>
 */

#include <sstream> 
#include "include/errors.hpp"

namespace Errors{
  void ConvergenceErr(std::ostringstream& output) {
    BaseErr("Convergence", output);
  }

  void ConvergenceErr(const char* output) {
    std::ostringstream pass_output;
    pass_output << output;
    BaseErr("Convergence", pass_output);
  }


}
