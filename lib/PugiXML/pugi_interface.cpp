#include "include/pugi_interface.h"
#include "Communicator/comm_io.hpp"
#include <string.h>
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <iterator>

namespace XML 
{
  

  int load_file(pugi::xml_document &doc, const char* filename) {
    pugi::xml_parse_result res = doc.load_file(filename);
    
    if (!res) {
      std::cout << "Error parsing XML input file ["<< filename <<"]\n";
      return 1;
    }
    return 0;
  }


  node select_node(pugi::xml_document &doc, const char *name) {
    return doc.child(name);
  }

  node getInputXML(const char* filename) {
    ParamsXML::instance();
    XML::node top_node; 
    
    XML::load_file(XML::ParamsXML::instance(), filename);
    //<Parameters> is a mandatory node for each input file
    top_node = XML::select_node(XML::ParamsXML::instance(),"Parameters");  
    return top_node;
  }

  bool attribute_compare(pugi::xml_node &node, 
			 const char* attribute_name,
			 const char* string){
    return strcmp(node.attribute(attribute_name).value(),string);
  }

  void descend(pugi::xml_node &node) {
    node = node.first_child();
    if (node ==NULL) {
      std::cout << "No nodes to descend" << std::endl;
    }
  }

  void descend(pugi::xml_node &node, const char *name) {
    node = node.child(name);
  }
  
  void next_sibling(pugi::xml_node &node, const char *name){
    node = node.next_sibling(name);
  }


  void read(pugi::xml_node node, const char *name, int& value, bool type) {
    if (node.child(name)!=NULL){
    value = atoi(node.child_value(name));
    } else {
      if (type == MANDATORY) {
	CCIO::cout << "Error: mandatory node ["<< name << "] not found\n"; 
	abort();
      } else {
	CCIO::cout << "Warning: node ["<< name << "] not found\n"; 
      }
    }
  }

  void read(pugi::xml_node node, const char *name, unsigned long& value, bool type) {
    char *end;
    value = strtol(node.child_value(name),&end, 0);

  }


  void read(pugi::xml_node node, const char *name, double& value) {
    value = atof(node.child_value(name));
  }

  void read(pugi::xml_node node, const char *name, bool& value) {
 
  }


  void read_array(pugi::xml_node node, 
		  const char* name, 
		  std::vector<int>& array) {
    //array is assumed to be of the correct size
    using namespace std;
    vector<string> tokens;
    int it;
    string read_res;
    
    
    read_res = node.child_value(name); 
    istringstream iss(read_res);
    copy(istream_iterator<string>(iss),
	 istream_iterator<string>(),
	 back_inserter<vector<string> >(tokens));

   if (array.size()){// only if array size was set before        
     if (array.size()!=tokens.size()) {
       cerr << "Check [" << name << "] number of parameters. Expected : " 
	    << array.size() << " found : " <<tokens.size() <<  endl;
       abort();}
   } else {
     array.resize(tokens.size());
   }
    for (it = 0; it<tokens.size(); it++) {
      array[it] = atoi(tokens[it].c_str());

    }


  }

  void read_array(pugi::xml_node node, 
		  const char* name, 
		  std::vector<unsigned long>& array) {
    //array is assumed to be of the correct size
    using namespace std;
    vector<string> tokens;
    int it;
    string read_res;
    char *end;
    
    read_res = node.child_value(name); 
    istringstream iss(read_res);
    copy(istream_iterator<string>(iss),
	 istream_iterator<string>(),
	 back_inserter<vector<string> >(tokens));

   if (array.size()){// only if array size was set before    
    if (array.size()!=tokens.size()) {
      cerr << "Check [" << name << "] number of parameters. Expected : " 
	   << array.size() << " found : " <<tokens.size() <<  endl;
      abort();}
   }else {
     array.resize(tokens.size());
   }

    for (it = 0; it<tokens.size(); it++) {
      array[it] = strtol(tokens[it].c_str(),&end,0);

    }
   

  }


  void read_array(pugi::xml_node node, 
		  const char* name, 
		  std::vector<double>& array) {
    //array is assumed to be of the correct size
    using namespace std;
    vector<string> tokens;
    int it;
    string read_res;
    
    
    read_res = node.child_value(name); 
    istringstream iss(read_res);
    copy(istream_iterator<string>(iss),
	 istream_iterator<string>(),
	 back_inserter<vector<string> >(tokens));
  
    if (array.size()){// only if array size was set before
      if (array.size()!=tokens.size()) {
	cerr << "Check [" << name << "] number of parameters. Expected : " 
	     << array.size() << " found : " <<tokens.size() <<  endl;
	abort();}
    } else {
     array.resize(tokens.size());
   }
      
      for (it = 0; it<tokens.size(); it++) {
	array[it] = atof(tokens[it].c_str());
	
      }
    

  }



};

