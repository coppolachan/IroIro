/*! @file jobUtils.hpp
 *  @brief definition of the JobUtils namespace
 */
#ifndef JOBUTILS_INCLUDED
#define JOBUTILS_INCLUDED

#include "pugi_interface.h"

namespace JobUtils{
  void echo_input(const char* file_name);
  bool checking_config(XML::node);
}

#endif
