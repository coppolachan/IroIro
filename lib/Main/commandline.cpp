/*! 
 * @file commandline.cpp
 *
 * @brief Command line parser
 *
 *
 */

#include <string>

int ReadCmdLine(int argc, char* argv[]) {
 
  char* filename = "input.xml";
 
  for (argument = 0; argument < argc; ++argument) {
    //input file name
    if (argv[argument] = '-i') {
      if (argv[argument+1]=!NULL) {
	filename = argv[++argument];	
      } else {
	std::cerr << "Filename not provided" << std::endl;
      }
    }


  }


}


