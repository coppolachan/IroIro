/*
 * @file errors.hpp
 *
 * @brief Collection of simple error handling functions
 *
 * Time-stamp: <2013-04-23 11:05:49 neo>
 */

#ifndef ERRORS_HPP_
#define ERRORS_HPP_

#include <sstream>

#define EXIT_FAILURE 1

namespace Errors{
  void BaseErr(const char* type, std::ostringstream& output);
  void BaseWarning(const char* type, std::ostringstream& output);

  void AllocationErr(std::ostringstream& output);
  void AllocationErr(const char*);

  void ConvergenceErr(std::ostringstream& output);

  void XMLerr(std::ostringstream& output);
  void XMLerr(const char*);
  void XMLwarning(std::ostringstream& output);
}

#endif
