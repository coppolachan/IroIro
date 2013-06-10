/*! 
 * @file ILDG_reader.hpp 
 *
 * @brief Declaration of the ILDGReader 
 *
 * Time-stamp: <2013-06-04 17:23:49 neo>
 */
#ifndef BINARY_READER_HPP_
#define BINARY_READER_HPP_

#include "general_reader.hpp"

namespace CCIO {
  class ILDGReader: public GeneralReader {
    // no header
    FILE* inputFile;
    unsigned int total_volume;
    unsigned int total_internal;
    unsigned int total_external;
  public:
    void set_sources(FILE *);
    int format(int gsite, int in, int ex) const;
    int read(double *buffer, unsigned int size);
    // QCDheader contains infos about the
    // source of data (FILE or LimeReader for example)
    int header();//reads the header, if present
    ILDGReader(int tot_vol,int tot_in,int tot_ex);
  };

}


#endif
