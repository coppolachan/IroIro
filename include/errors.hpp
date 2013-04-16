/*
 * @file errors.hpp
 *
 * @brief Collection of simple error handling functions
 *
 * Time-stamp: <2013-04-16 14:35:30 neo>
 */

#ifndef ERRORS_HPP_
#define ERRORS_HPP_

#include <sstream>

#define EXIT_FAILURE 1

namespace Errors{
  void XMLerr(std::ostringstream& output);
  void XMLwarning(std::ostringstream& output);
}

#endif
