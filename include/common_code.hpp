/*!
 * @file common_code.hpp
 *
 * @brief Common code general include file
 *
 * Include this file into your main code.
 *
 * Includes the fundamental features common to every code:
 * - Fields and common parameters
 * - Communicator headers
 * - XML interface
 */

#ifndef COMMON_CODE_HPP_
#define COMMON_CODE_HPP_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

//Fields and common parameters
#include "include/common_fields.hpp"
#include "include/commonPrms.h"

//Communicator headers
#include "Communicator/communicator.h"
#include "Communicator/comm_io.hpp"
#include "Communicator/fields_io.hpp"

//XML interface
#include "include/pugi_interface.h"

#endif


