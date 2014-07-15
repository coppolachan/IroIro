/*! 
 * @file text_header.hpp 
 *
 * @brief Declaration of the Text_header class
 *
 * Time-stamp: <2014-07-15 12:56:41 neo>
 */
#ifndef TEXT_HEADER_HPP_
#define TEXT_HEADER_HPP_

#include "generic_header.hpp"
#include <string>
#include <stdint.h>

namespace CCIO {

 class Text_header: public QCDheader {
    std::map<std::string, std::string> tokens;
  public:
   void add_element(std::string, std::string);
   void get_value(std::string, int&);
   void get_value(std::string, uint32_t&);
   void get_value(std::string, uint64_t&);
   void get_value(std::string, double& );
   void get_value(std::string, std::string&);
   void print();
   Text_header();

  };


}


#endif
