/*! 
 * @file NERSC_writer.hpp 
 *
 * @brief Declaration of the NERSCWriter 
 *
 * Time-stamp: <2014-07-14 17:21:29 neo>
 */
#ifndef NERSC_WRITER_HPP_
#define NERSC_WRITER_HPP_

#include "general_writer.hpp"
#include "include/field.h"

namespace CCIO {
  class NERSCWriter: public GeneralWriter {
    FILE* outputFile;
    unsigned int num_blocks;
    unsigned int total_volume;
    unsigned int total_internal;
    unsigned int total_external;
    unsigned int writer_first_call;
    double plaquette;
  public:
    void set_output(FILE *);
    int format(int gsite, int in, int ex) const;
    int write(double *buffer, unsigned int size);
    int header();
    NERSCWriter(int tot_vol,int tot_in,int tot_ex, const Field& U);
  };

}


#endif
