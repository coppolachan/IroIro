/*! 
 * @file ILDG_reader.cpp 
 *
 * @brief Definition of the ILDGReader class methods
 *
 * Time-stamp: <2013-06-06 12:00:43 neo>
 */

#include "ILDG_reader.hpp"
#include "include/messages_macros.hpp"
#include "Tools/byteswap.hpp"
#include <string.h>

namespace CCIO {
  void ILDGReader::set_sources(FILE *src){
    inputFile = src;
  }
  int ILDGReader::format(int gsite, int in, int ex) const{
      return in +total_internal*(ex+ total_external*gsite);
  }

  int ILDGReader::read(double *buffer, unsigned int size){
    _Message(DEBUG_VERB_LEVEL, "ILDG Reader...\n");
    //Read just as it is
    //order is (in, ex, sites)  
    n_uint64_t bytes_to_read = sizeof(float)*size;
    n_uint64_t read =  bytes_to_read;
    float uin[100000];
    limeReaderReadData(uin, &read, LimeR);

#ifndef BIG_ENDIAN_TYPE
    byte_swap(uin, size);
#endif


#ifdef DEBUG_VERB_LEVEL
    for (int i = 0; i < bytes_to_read; i++){
      std::cout << "buffer["<<i<<"] = "<< buffer[i] << "\n";
    }
#endif 
  }

  int ILDGReader::header(){
    _Message(DEBUG_VERB_LEVEL, "ILDG Header...\n");
     n_uint64_t nbytes;

    LimeR = limeCreateReader(inputFile);
    do {limeReaderNextRecord(LimeR);}
    while (strncmp(limeReaderType(LimeR), "ildg-binary-data",16));

    nbytes = limeReaderBytes(LimeR);//size of this record (configuration)
    _Message(DEBUG_VERB_LEVEL, "ILDG Header: reading "<<nbytes<< " bytes record...\n");
  }

  int ILDGReader::check(Field& U){};

  ILDGReader::ILDGReader(int tot_vol,int tot_in,int tot_ex):
    total_volume(tot_vol),total_internal(tot_in),
    total_external(tot_ex){};
}


