/*! 
 * @file commandline.hpp
 *
 * @brief Command line parser
 *
 *
 */

struct CommandOptions{
  char* filename;
  char* output;

};

CommandOptions ReadCmdLine(int argc, char* argv[]);
