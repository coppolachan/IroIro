/*
 * @file io_errors.cpp
 *
 * @brief Collection of simple I/O error handling functions
 *
 * Time-stamp: <2013-04-23 17:26:08 neo>
 */

#include <sstream> 
#include "include/errors.hpp"

namespace Errors{
  void IOErr(std::ostringstream& output) {
    BaseErr("I/O", output);
  }
  
  void IOErr(IOtype type, const char* MoreInfo = 0) {
    std::ostringstream msg;
    switch(type) {
    case FileNotFound:
      msg << "File not found: " << MoreInfo;
      break;
    default: 
      msg << "Generic I/O error.";
      break;
    }
    IOErr(msg);
  }
 
 

}
