/*! 
 * @file binary_writer.cpp 
 *
 * @brief Definition of the BinaryWriter class methods
 *
 * Time-stamp: <2013-09-12 15:18:44 cossu>
 */

#include "binary_writer.hpp"

namespace CCIO {
  void BinaryWriter::set_output(FILE *out){
    outputFile = out;
  }
  int BinaryWriter::format(int gsite, int in, int ex) const{
      return in +total_internal*(ex+ total_external*gsite);
  }

  int BinaryWriter::write(double *buffer, unsigned int size){
    size_t res = fwrite(buffer, sizeof(double) ,size , outputFile);
    return res;
  }

  int BinaryWriter::header(){};

  BinaryWriter::BinaryWriter(int tot_vol,int tot_in,int tot_ex):
    total_volume(tot_vol),total_internal(tot_in),
    total_external(tot_ex){};
}


