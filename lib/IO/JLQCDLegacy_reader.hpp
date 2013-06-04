/*! 
 * @file JLQCDLegacy_reader.hpp 
 *
 * @brief Declaration of the JLQCDLegacyReader 
 *
 * Time-stamp: <2013-06-04 17:44:37 neo>
 */
#ifndef JLQCD_READER_HPP_
#define JLQCD_READER_HPP_

#include "general_reader.hpp"

namespace CCIO {
  class JLQCDLegacyReader: public GeneralReader {
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
    JLQCDLegacyReader(int tot_vol,int tot_in,int tot_ex);
  };

}


#endif
