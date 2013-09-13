/*! 
 * @file general_writer.hpp 
 *
 * @brief Declaration of the GeneralWriter abstract class
 *
 * Time-stamp: <2013-09-12 17:10:37 cossu>
 */
#ifndef GENERAL_WRITER_HPP_
#define GENERAL_WRITER_HPP_

#include <stdio.h>

namespace CCIO {
  class GeneralWriter {
    //global info at contruction time
  public:
    virtual void set_output(FILE*) = 0;
    virtual int format(int gsite, int in, int ex) const = 0; 
    virtual int write(double *buffer, unsigned int size) = 0;
    // QCDheader contains infos about the
    // source of data (FILE or LimeReader for example)
    virtual int header() = 0;//writes the header, if present
  };

}


#endif

