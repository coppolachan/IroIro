/*! 
 * @file binary_reader.hpp 
 *
 * @brief Declaration of the BinaryReader 
 *
 * Time-stamp: <2013-06-04 17:45:21 neo>
 */
#ifndef BINARY_READER_HPP_
#define BINARY_READER_HPP_

#include "general_reader.hpp"

namespace CCIO {
  class BinaryReader: public GeneralReader {
    // no header
    FILE* inputFile;
    unsigned int total_volume;
    unsigned int total_internal;
    unsigned int total_external;
  public:
    void set_sources(FILE *);
    int format(int gsite, int in, int ex) const;
    int read(double *buffer, unsigned int size);
    int header();
    int check(Field&);
    BinaryReader(int tot_vol,int tot_in,int tot_ex);
  };

}


#endif
