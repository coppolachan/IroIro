/*! 
 * @file generic_header.hpp 
 *
 * @brief Declaration of the QCDheader abstract class
 *
 * Time-stamp: <2013-06-05 11:33:57 neo>
 */
#ifndef GENERIC_HEADER_HPP_
#define GENERIC_HEADER_HPP_

#include "include/config.h"
#include <stdio.h>
#include <string>
#include <map>

namespace CCIO {
  //Abstract class for headers
  class QCDheader {
  public:
    virtual void get_value(std::string, int& ) = 0;
    virtual void get_value(std::string, double& ) = 0;
    virtual void get_value(std::string, std::string& ) = 0;
    virtual void print() = 0;
  };


} //end of namespace

#endif
