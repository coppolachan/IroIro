/*
 * @file allocation_errors.cpp
 *
 * @brief Collection of simple allocation error handling functions
 *
 * Time-stamp: <2013-10-23 16:48:12 noaki>
 */

#include <stdlib.h> // exit definition
#include "include/errors.hpp"
#include "Communicator/comm_io.hpp"

namespace Errors{
  void BaseErr(const char* type, std::ostringstream& output) {
    CCIO::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    CCIO::cerr << "["<<type<<" ERROR] Description:\n"<< output.str() << "\n";
    CCIO::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    CCIO::cerr << "Exiting...\n";
    exit(EXIT_FAIL);
  }

  void BaseWarning(const char* type, std::ostringstream& output) {
    CCIO::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    CCIO::cerr << "["<<type<<" WARNING] Description:\n"<< output.str() << "\n";
    CCIO::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  }


}
