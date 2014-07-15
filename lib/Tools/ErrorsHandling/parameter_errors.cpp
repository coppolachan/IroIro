/*
 * @file parameter_errors.cpp
 *
 * @brief Collection of parameter errors handling functions
 *
 * Time-stamp: <2014-07-15 12:09:51 neo>
 */

#include <sstream> 
#include "include/errors.hpp"

namespace Errors{
 
  void ParameterErr(std::ostringstream& output) {
    BaseErr("Parameter", output);
  }

  void ParameterErr(const char* message) {
    std::ostringstream pass_msg;
    pass_msg << message;
    BaseErr("Parameter", pass_msg);
  }

  void ParameterWarning(std::ostringstream& output) {
    BaseWarning("Parameter", output);
  }


}
