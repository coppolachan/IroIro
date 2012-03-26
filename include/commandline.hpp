/*! 
 * @file commandline.hpp
 * @brief Command line parser
 */

#ifndef COMMANDLINE_H_
#define COMMANDLINE_H_

struct CommandOptions{
  char* filename;
  char* output;
};

CommandOptions ReadCmdLine(int argc, char* argv[]);

#endif
