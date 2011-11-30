#ifndef PUGI_INTERFACE_H_
#define PUGI_INTERFACE_H_

#include "PugiXML/pugixml.hpp"
#include "include/singleton.h"
#include <vector>

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

  void descend(pugi::xml_node &node, const char *name);

  void next_sibling(pugi::xml_node &node, const char *name);

  bool attribute_compare(pugi::xml_node &node, 
			 const char* attrname, 
			 const char* string);

  void read(pugi::xml_node node, const char *name, int& value, bool type = false);

  void read(pugi::xml_node node, const char *name, unsigned long& value, bool type = false);

  void read(pugi::xml_node node, const char *name, double& value);

  void read(pugi::xml_node node, const char *name, bool& value);

  void read_array(pugi::xml_node node, 
		  const char* name, 
		  std::vector<int>& array);

  void read_array(pugi::xml_node node, 
		  const char* name, 
		  std::vector<unsigned long>& array);


  void read_array(pugi::xml_node node, 
		  const char* name, 
		  std::vector<double>& array);




};


#endif//PUGI_INTERFACE_H_
