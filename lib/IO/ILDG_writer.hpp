/*! 
 * @file ILDG_writer.hpp 
 *
 * @brief Declaration of the ILDGWriter 
 *
 * Time-stamp: <2013-09-12 17:11:00 cossu>
 */
#ifndef ILDG_WRITER_HPP_
#define ILDG_WRITER_HPP_

#include "general_writer.hpp"
extern "C" { // for linkage
#include "lime.h"
}

namespace CCIO {
  class ILDGWriter: public GeneralWriter {
    FILE* outputFile;
    LimeWriter* LimeW;
    LimeRecordHeader* LimeHeader;
    unsigned int num_blocks;
    unsigned int total_volume;
    unsigned int total_internal;
    unsigned int total_external;
    unsigned int writer_first_call;
  public:
    void set_output(FILE *);
    int format(int gsite, int in, int ex) const;
    int write(double *buffer, unsigned int size);
    int header();
    ILDGWriter(int tot_vol,int tot_in,int tot_ex);
  };

}


#endif
