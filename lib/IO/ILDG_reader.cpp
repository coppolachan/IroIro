/*! 
 * @file ILDG_reader.cpp 
 *
 * @brief Definition of the ILDGReader class methods
 *
 * Time-stamp: <2013-06-04 17:27:14 neo>
 */


#include "binary_reader.hpp"
#include "lime.h"

namespace CCIO {
  void ILDGReader::set_sources(FILE *src){
    inputFile = src;
  }
  int ILDGReader::format(int gsite, int in, int ex) const{
      return in +total_internal*(ex+ total_external*gsite);
  }

  int ILDGReader::read(double *buffer, unsigned int size){
    //Read just as it is
    //order is (in, ex, sites)  
    size_t res = fread(buffer,sizeof(double),size,inputFile);
    return res;
  }

  int ILDGReader::header(){
    LimeReader *reader;
    n_uint64_t nbytes;

    reader = limeCreateReader(inputFile);
    do {limeReaderNextRecord(reader);}
    while (strncmp(limeReaderType(reader), "ildg-binary-data",16));

    nbytes = limeReaderBytes(reader);//size of this record (configuration)
  }

  };

  ILDGReader::ILDGReader(int tot_vol,int tot_in,int tot_ex):
    total_volume(tot_vol),total_internal(tot_in),
    total_external(tot_ex){};
}


