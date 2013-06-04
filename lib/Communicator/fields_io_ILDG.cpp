/*! 
 * @file fields_io_ILDG.cpp 
 *
 * @brief Declarations of MPI safe ILDG format reading routines for fields
 *
 * Time-stamp: <2013-06-04 10:49:40 neo>
 */

#include "fields_io.hpp"

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "lime.h"


namespace CCIO {
  QCDheader* ReadILDGheader(FILE *inputFile, fpos_t *pos){
    LimeReader *reader;
    n_uint64_t nbytes;

    reader = limeCreateReader(inputFile);
    do {limeReaderNextRecord(reader);}
    while (strncmp(limeReaderType(reader), "ildg-binary-data",16));

    nbytes = limeReaderBytes(reader);//size of this record (configuration)
  }


}
