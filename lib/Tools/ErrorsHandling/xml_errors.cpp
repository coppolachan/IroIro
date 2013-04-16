/*
 * @file xml_errors.cpp
 *
 * @brief Collection of simple XML error handling functions
 *
 * Time-stamp: <2013-04-16 14:35:41 neo>
 */

#include <stdlib.h> // exit definition
#include "include/errors.hpp"
#include "Communicator/comm_io.hpp"

namespace Errors{
 
  void XMLerr(std::ostringstream& output) {
    CCIO::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    CCIO::cerr << "[XML ERROR] Description:\n"<< output.str() << "\n";
    CCIO::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    CCIO::cerr << "Exiting...\n";
    exit(EXIT_FAILURE);
  }

  void XMLwarning(std::ostringstream& output) {
    CCIO::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    CCIO::cout << "[XML Warning] Description:\n"<< output.str() << "\n";
    CCIO::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  }


}
