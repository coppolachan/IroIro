/*! 
 * @file commandline.hpp
 * @brief Command line parser declaration
 *
 * Time-stamp: <2013-04-25 17:01:16 neo>
 */

#ifndef COMMANDLINE_H_
#define COMMANDLINE_H_

struct CommandOptions{
  char* filename;
  char* output;
};

CommandOptions ReadCmdLine(int argc, char* argv[]);

#endif
