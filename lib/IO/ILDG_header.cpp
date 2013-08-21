/*! 
 * @file ILDG_header.cpp 
 *
 * @brief Definition of the ILDG_header class
 *
 * Time-stamp: <2013-06-24 15:16:26 neo>
 */

#include "ILDG_header.hpp"

namespace CCIO {
   
  void ILDG_header::fill_xml(void* xmlbuf, size_t bytes){
    XML::load_buffer(ILDG_XML,xmlbuf, bytes);
    ILDGnode = XML::select_node(ILDG_XML, "ildgFormat");
      /*
    XML::read(ILDGnode, "precision", ILDG_h.precision, MANDATORY);
    XML::read(ILDGnode, "version", ILDG_h.version, MANDATORY);
    XML::read(ILDGnode, "lx", ILDG_h.lx, MANDATORY);
    XML::read(ILDGnode, "ly", ILDG_h.ly, MANDATORY);
    XML::read(ILDGnode, "lz", ILDG_h.lz, MANDATORY);
    XML::read(ILDGnode, "lt", ILDG_h.lt, MANDATORY);
    */

    
  };

  void ILDG_header::add_element(std::string, std::string){};  
  void ILDG_header::get_value(std::string name, int& value){
    XML::read(ILDGnode, name.c_str(), value, MANDATORY);
  };
  void ILDG_header::get_value(std::string name, double& value){
    XML::read(ILDGnode, name.c_str(), value, MANDATORY);
  };
  void ILDG_header::get_value(std::string name, std::string& value){
    XML::read(ILDGnode, name.c_str(), value, MANDATORY);
  };
  void ILDG_header::print(){};
  


}
