/*! 
 * @file ILDG_reader.hpp 
 *
 * @brief Declaration of the ILDGReader 
 *
 * Time-stamp: <2013-06-05 10:41:56 neo>
 */
#ifndef ILDG_READER_HPP_
#define ILDG_READER_HPP_

#include "general_reader.hpp"
extern "C" { // for linkage
#include "lime.h"
}

namespace CCIO {
  class ILDGReader: public GeneralReader {
    FILE* inputFile;
    LimeReader* LimeR;
    unsigned int total_volume;
    unsigned int total_internal;
    unsigned int total_external;
  public:
    void set_sources(FILE *);
    int format(int gsite, int in, int ex) const;
    int read(double *buffer, unsigned int size);
    int header();
    int check(Field&);
    ILDGReader(int tot_vol,int tot_in,int tot_ex);
  };

}


#endif
