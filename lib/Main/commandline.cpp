/*! 
 * @file commandline.cpp
 *
 * @brief Command line parser
 *
 * Time-stamp: <2013-04-26 16:54:05 neo>
 */
#include "include/commandline.hpp"
#include <string>
#include <iostream>
#include <algorithm>

char* getCmdOption(char ** begin, char ** end, const std::string & option){
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end) return *itr;
  return 0;
}

CommandOptions ReadCmdLine(int argc, char* argv[]) {
  CommandOptions Result;
  Result.filename = getCmdOption(argv, argv + argc, "-i");

  //default input XML
  if(!Result.filename) Result.filename = const_cast<char*>("input.xml"); 
  
  Result.output = const_cast<char*>("output.xml");
  return Result;
}

