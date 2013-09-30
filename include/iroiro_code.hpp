/*!
 * @file common_code.hpp
 *
 * @brief IroIro code general include file
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

#include "include/iroiro_config.h"

//Fields and common parameters
#include "include/common_fields.hpp"
#include "include/commonPrms.hpp"
#include "EigenModes/eigenModes.hpp"

//Communicator headers
#include "Communicator/communicator.hpp"
#include "Communicator/comm_io.hpp"

//IO headers
#include "IO/fields_io.hpp"

//XML interface
#include "include/pugi_interface.h"

//Numerical constants
#include "include/numerical_const.hpp"




#endif


