/*
 * @file allocation_errors.cpp
 *
 * @brief Collection of simple allocation error handling functions
 *
 * Time-stamp: <2013-04-22 13:09:12 cossu>
 */

#include <sstream> 
#include "include/errors.hpp"

namespace Errors{
  void AllocationErr(std::ostringstream& output) {
    BaseErr("Allocation", output);
  }

  void AllocationErr(const char* output) {
    std::ostringstream pass_output;
    pass_output << output;
    BaseErr("Allocation", pass_output);
  }


}
