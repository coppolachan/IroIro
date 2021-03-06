/*
 * @file errors.hpp
 *
 * @brief Collection of simple error handling functions
 *
 * Time-stamp: <2013-09-30 21:57:08 noaki>
 */

#ifndef ERRORS_HPP_
#define ERRORS_HPP_

#include <sstream>

#define EXIT_FAIL 1

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

  void ParameterErr(const char*);
  void ParameterErr(std::ostringstream& output);
  void ParameterWarning(std::ostringstream& output);

}

#endif
