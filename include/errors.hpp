/*
 * @file errors.hpp
 *
 * @brief Collection of simple error handling functions
 *
 * Time-stamp: <2013-06-03 14:27:03 neo>
 */

#ifndef ERRORS_HPP_
#define ERRORS_HPP_

#include <sstream>

#define EXIT_FAILURE 1

typedef std::ostringstream ErrorString;

namespace Errors{
  typedef enum IOerrorType {GenericError, FileNotFound, ParsingError} IOtype;

  void BaseErr(const char* type, std::ostringstream& output);
  void BaseWarning(const char* type, std::ostringstream& output);

  void AllocationErr(std::ostringstream& output);
  void AllocationErr(const char*);

  void ConvergenceErr(std::ostringstream& output);
  
  void IOErr(std::ostringstream& output);
  void IOErr(IOtype, const char*);

  void XMLerr(std::ostringstream& output);
  void XMLerr(const char*);
  void XMLwarning(std::ostringstream& output);
}

#endif
