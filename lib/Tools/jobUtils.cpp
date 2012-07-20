/*! @file jobUtils.cpp
 * @brief implementation of member functions of JobUtils
 */
#include "jobUtils.hpp"
#include "lib/Communicator/comm_io.hpp"
#include <fstream>

namespace JobUtils{
  void echo_input(const char* file_name){

    std::ifstream input(file_name);
    std::string buf;
    CCIO::cout<<"==== [begin] Contents of the input file "
	      << file_name << "===="
	      << std::endl;

    while(input && getline(input,buf))
      CCIO::cout<< buf << std::endl;

    CCIO::cout<<"==== [end] Contents of the input file "
	      << file_name << "===="
	      << std::endl;
    return;
  }
}
