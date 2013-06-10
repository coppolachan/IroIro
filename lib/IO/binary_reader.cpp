/*! 
 * @file binary_reader.cpp 
 *
 * @brief Definition of the BinaryReader class methods
 *
 * Time-stamp: <2013-06-04 17:49:54 neo>
 */


#include "binary_reader.hpp"


namespace CCIO {
  void BinaryReader::set_sources(FILE *src){
    inputFile = src;
  }
  int BinaryReader::format(int gsite, int in, int ex) const{
      return in +total_internal*(ex+ total_external*gsite);
  }

  int BinaryReader::read(double *buffer, unsigned int size){
    //Read just as it is
    //order is (in, ex, sites)  
    size_t res = fread(buffer,sizeof(double),size,inputFile);
    return res;
  }

  int BinaryReader::header(){};
  int BinaryReader::check(Field&){}; // no checks

  BinaryReader::BinaryReader(int tot_vol,int tot_in,int tot_ex):
    total_volume(tot_vol),total_internal(tot_in),
    total_external(tot_ex){};
}


