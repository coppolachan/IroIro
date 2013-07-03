/*! 
 * @file general_reader.hpp 
 *
 * @brief Declaration of the GeneralReader abstract class
 *
 * Time-stamp: <2013-06-04 17:44:16 neo>
 */
#ifndef GENERAL_READER_HPP_
#define GENERAL_READER_HPP_

#include <stdio.h>
#include "include/field.h"

namespace CCIO {
  class GeneralReader {
    //global info at contruction time
  public:
    virtual void set_sources(FILE*) = 0;
    virtual int format(int gsite, int in, int ex) const = 0; 
    virtual int read(double *buffer, unsigned int size) = 0;
    // QCDheader contains infos about the
    // source of data (FILE or LimeReader for example)
    virtual int header() = 0;//reads the header, if present
    virtual int check(Field&) = 0; 
  };

}


#endif

