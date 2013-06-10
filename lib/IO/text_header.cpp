/*! 
 * @file text_header.cpp 
 *
 * @brief Definition of the Text_header class
 *
 * Time-stamp: <2013-06-04 17:09:54 neo>
 */

#include "text_header.hpp"
#include <stdlib.h>
#include <iostream>

namespace CCIO {
  Text_header::Text_header(){
  };

  void Text_header::add_element(std::string key, std::string value){
    tokens.insert(std::pair<std::string, std::string>(key, value));
  }

  void Text_header::get_value(std::string key, int& value){
    value = atoi(tokens[key].c_str());
  }

  void Text_header::get_value(std::string key, double& value) {
    value = atof(tokens[key].c_str());
  }

  void Text_header::get_value(std::string key, std::string& value){
    value.assign(tokens[key].c_str());
  }

  void Text_header::print(){
    std::map<std::string, std::string>::iterator it;
    for (it=tokens.begin(); it!=tokens.end(); ++it)
    std::cout << it->first << " => " << it->second << '\n';
  }

}
