/*! 
 * @file NERSC_reader.hpp 
 *
 * @brief Declaration of the NERSCReader 
 *
 * Time-stamp: <2014-07-15 11:29:20 neo>
 */
#ifndef NERSC_READER_HPP_
#define NERSC_READER_HPP_

#include "general_reader.hpp"
#include "text_header.hpp"

namespace CCIO {
  class NERSCReader: public GeneralReader {
    // no header
    FILE* inputFile;
    Text_header TH;
    unsigned int total_volume;
    unsigned int total_internal;
    unsigned int total_external;
    unsigned int bits;
    bool is3x2;
    double stored_plaquette;
    double stored_link;
    double tolerance;

    void NERSCtoIROIRO_f(double* out, int block);
    void NERSCtoIROIRO_d(double* out, int block);
    template <typename FLOAT>
    void reconstruct3x3(double* out, FLOAT* in,int block3x2);
  public:
    void set_sources(FILE *);
    int format(int gsite, int in, int ex) const;
    int read(double *buffer, unsigned int size);
    // QCDheader contains infos about the
    // source of data (FILE or LimeReader for example)
    int header();//reads the header, if present
    int check(Field&);
    NERSCReader(int tot_vol,int tot_in,int tot_ex);
  };

}



#endif
