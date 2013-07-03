/*! 
 * @file ILDG_header.hpp 
 *
 * @brief Declaration of the ILDG_header class
 *
 * Time-stamp: <2013-06-04 15:42:20 neo>
 */
#ifndef ILDG_HEADER_HPP_
#define ILDG_HEADER_HPP_

#include "generic_header.hpp"
#include <string>

namespace CCIO {

 class ILDG_header: public QCDheader {
 
  public:
    void add_element(std::string, std::string);
    void get_value(std::string, int&);
    void get_value(std::string, double& );
    void get_value(std::string, std::string& );
    void print();
    ILDG_header();

  };


}


#endif
