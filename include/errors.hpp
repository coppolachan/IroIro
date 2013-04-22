/*
 * @file errors.hpp
 *
 * @brief Collection of simple error handling functions
 *
 * Time-stamp: <2013-04-22 13:06:14 cossu>
 */

#ifndef ERRORS_HPP_
#define ERRORS_HPP_

#include <sstream>

#define EXIT_FAILURE 1

namespace Errors{
  void BaseErr(const char* type, std::ostringstream& output);
  void AllocationErr(std::ostringstream& output);
  void AllocationErr(const char*);
  void XMLerr(std::ostringstream& output);
  void XMLwarning(std::ostringstream& output);
}

#endif
