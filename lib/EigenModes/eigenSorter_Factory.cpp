/*! @file eigenSorter_Factory.cpp
 *  @brief eigenSorter-factories selector 
 */
#include "eigenSorter_Factory.hpp"
#include "Communicator/comm_io.hpp"
#include <string.h>
#include <cstdlib>

namespace EigenModes{
  EigenSorterFactory* createEigenSorterFactory(XML::node node){
    if(node != NULL){
      std::string sort_type;
      XML::read(node,"Sorter",sort_type,MANDATORY);

      if(sort_type =="LowestModes")  return new EigenSorterFactory_low();
      if(sort_type =="HighestModes") return new EigenSorterFactory_high();

      CCIO::cerr<<"No EigenSorter available with name["
		<< sort_type << "]" << std::endl;
      abort();
    }else{
      CCIO::cout << "Requested node is missing in input file "
                 << "(EigenSorter Object)\n"
                 << "Request by " << node.parent().name() << std::endl;
      abort();
    }
  }
}
