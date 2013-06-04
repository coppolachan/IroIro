/*! 
 * @file JLQCD_reader.cpp 
 *
 * @brief Definition of the JLQCDReader class methods
 *
 * Time-stamp: <2013-06-04 17:49:33 neo>
 */


#include "JLQCDLegacy_reader.hpp"


namespace CCIO {
  void JLQCDLegacyReader::set_sources(FILE *src){
    inputFile = src;
  }
  int JLQCDLegacyReader::format(int gsite, int in, int ex) const{
    return in +total_internal*(gsite+ total_volume*ex);
  }

  int JLQCDLegacyReader::read(double *buffer, unsigned int size){
    //order is (in, sites, ex)
    //Read #tot_ex blocks
    int chunk = size/total_external;
    fpos_t pos;
    fgetpos(inputFile, &pos);
    for(int ext=0; ext<total_external; ++ext){
      //reads the ex blocks
      size_t res = fread((buffer+ext*chunk), 
			 sizeof(double), 
			 chunk, 
			 inputFile);
      fseek(inputFile, sizeof(double)*total_volume*total_internal , SEEK_CUR);
    }
    fsetpos(inputFile, &pos);
    fseek(inputFile, sizeof(double)*chunk, SEEK_CUR);

  }

  int JLQCDLegacyReader::header(){};//no header
  int JLQCDLegacyReader::check(Field&){}; // no checks

  // Constructor
  JLQCDLegacyReader::JLQCDLegacyReader(int tot_vol,int tot_in,int tot_ex):
    total_volume(tot_vol),total_internal(tot_in),
    total_external(tot_ex){};
}


