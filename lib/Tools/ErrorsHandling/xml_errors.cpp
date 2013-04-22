/*
 * @file xml_errors.cpp
 *
 * @brief Collection of simple XML error handling functions
 *
 * Time-stamp: <2013-04-22 16:49:24 neo>
 */

#include <sstream> 
#include "include/errors.hpp"

namespace Errors{
 
  void XMLerr(std::ostringstream& output) {
    BaseErr("XML", output);
  }

  void XMLerr(const char* message) {
    std::ostringstream pass_msg;
    pass_msg << message;
      BaseErr("XML", pass_msg);
  }

  void XMLwarning(std::ostringstream& output) {
    BaseWarning("XML", output);
  }


}
