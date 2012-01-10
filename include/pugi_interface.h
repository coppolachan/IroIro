#ifndef PUGI_INTERFACE_H_
#define PUGI_INTERFACE_H_

#include "PugiXML/pugixml.hpp"
#include "include/singleton.h"
#include <vector>
#include <string>
#define MANDATORY true

/*!
 * @brief The namespace for the %XML interface functions
 *
 * Uses <a href="http://pugixml.org/">pugixml</a> parser the the core functions.
 */
namespace XML 
{
  typedef pugi::xml_document document;
  typedef pugi::xml_node node;

  typedef Singleton<document> ParamsXML;

  int load_file(pugi::xml_document &doc, const char *filename);

  node select_node(pugi::xml_document &doc, const char *name);

  node getInputXML(const char *filename);

  void descend(pugi::xml_node &node);

  void descend(pugi::xml_node &node, const char *name, bool type = false);

  void next_sibling(pugi::xml_node &node, const char *name, bool type = false);

  bool attribute_compare(pugi::xml_node &node, 
			 const char* attrname, 
			 const char* string);

  int read(pugi::xml_node node, const char *name, int& value, bool type = false);

  int read(pugi::xml_node node, const char *name, unsigned long& value, bool type = false);

  int read(pugi::xml_node node, const char *name, double& value,  bool type = false);

  int read(pugi::xml_node node, const char *name, bool& value, bool type = false);

  int read(pugi::xml_node node, const char *name, std::string& string, bool type = false);


  int read_array(pugi::xml_node node, 
		  const char* name, 
		  std::vector<int>& array, bool type = false);

  int read_array(pugi::xml_node node, 
		  const char* name, 
		  std::vector<unsigned long>& array, bool type = false);


  int read_array(pugi::xml_node node, 
		  const char* name, 
		  std::vector<double>& array, bool type = false);




};


#endif//PUGI_INTERFACE_H_
