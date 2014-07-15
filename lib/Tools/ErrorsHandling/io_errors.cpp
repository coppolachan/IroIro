/*
 * @file io_errors.cpp
 *
 * @brief Collection of simple I/O error handling functions
 *
 * Time-stamp: <2014-07-15 12:10:50 neo>
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
      msg << "Generic I/O error.  " << MoreInfo;
      break;
    }
    IOErr(msg);
  }
 
 

}
