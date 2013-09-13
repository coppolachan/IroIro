/*! 
 * @file binary_writer.hpp 
 *
 * @brief Declaration of the BinaryWriter class 
 *
 * Time-stamp: <2013-09-12 15:18:25 cossu>
 */
#ifndef BINARY_WRITER_HPP_
#define BINARY_WRITER_HPP_

#include "general_writer.hpp"

namespace CCIO {
  class BinaryWriter: public GeneralWriter {
    // no header
    FILE* outputFile;
    unsigned int total_volume;
    unsigned int total_internal;
    unsigned int total_external;
  public:
    void set_output(FILE *);
    int format(int gsite, int in, int ex) const;
    int write(double *buffer, unsigned int size);
    int header();
    BinaryWriter(int tot_vol,int tot_in,int tot_ex);
  };
}



#endif
