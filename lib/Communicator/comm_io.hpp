/*!
 * @file comm_io.hpp 
 * @brief Declaration of MPI safe input/output routines
 *
 * Overrides the standard cout and cin in an independent namespace
 */

#ifndef COMMCODE_IO_H_
#define COMMCODE_IO_H_

#include <istream>
#include <ostream>
#include <iomanip>
#include <string>
#include <vector>
#include "include/macros.hpp"

//------------------------------------------------------------------
/*! @brief StandardInputStream class for parallel 
 * version of standard input
 */

class StandardInputStream
{
public:
  //! Constructor
  StandardInputStream();

  //! destructor
  ~StandardInputStream();

  //! Constructor from a stream
  void init(std::istream *is);

//  //! Redirect input stream
//  void rdbuf(std::streambuf* b);

  //! Return true if some failure occurred in previous IO operation
  bool fail();

  // Overloaded input functions
  StandardInputStream& operator>>(std::string& input);
  StandardInputStream& operator>>(char& input);
  StandardInputStream& operator>>(int& input);
  StandardInputStream& operator>>(unsigned int& input);
  StandardInputStream& operator>>(short int& input);
  StandardInputStream& operator>>(unsigned short int& input);
  StandardInputStream& operator>>(long int& input);
  StandardInputStream& operator>>(unsigned long int& input);
  StandardInputStream& operator>>(long long int& input);
  StandardInputStream& operator>>(float& input);
  StandardInputStream& operator>>(double& input);
  StandardInputStream& operator>>(long double& input);
  StandardInputStream& operator>>(bool& input);


private:
  //! Hide copy constructor and =
  StandardInputStream(const StandardInputStream&) {}
  void operator=(const StandardInputStream&) {}

protected:
  // The universal data-reader. All the read functions call this
  template<typename T>
  StandardInputStream& readPrimitive(T& output);

  // Get the internal istream
  std::istream& getIstream() {return *is;}

  //! Is the stream open?
  bool is_open() {return open;}

private:
  std::istream* is;
  bool open;
};


//! Read an array
template<class T>
inline
StandardInputStream& operator>>(StandardInputStream& s, std::vector<T>& d)
{
  for(int i=0; i < d.size(); ++i)
    s >> d[i];

  return s;
}



//-------------------------------------------------------------------
//! StandardOutputStream class
/*! Parallel version of standard output
 */
class StandardOutputStream
{
public:
  //! Constructor
  StandardOutputStream();

  ~StandardOutputStream();

  //! Constructor from a stream
  void init(std::ostream *os);

//  //! Redirect output stream
//  void rdbuf(std::streambuf* b);

  //! Flush the buffer
  void flush();

  //! Return true if some failure occurred in previous IO operation
  bool fail();

  //! Hook for manipulators
  StandardOutputStream& operator<<(std::ostream& (*op)(std::ostream&));
  StandardOutputStream& operator<<(StandardOutputStream& (*op)(StandardOutputStream&));
#ifndef HITACHISR16K
  friend StandardOutputStream& operator<<(StandardOutputStream&, std::_Setw f);
  friend StandardOutputStream& operator<<(StandardOutputStream&, std::_Setprecision f);
#endif
#ifdef HITACHISR16K //Hitachi hack
  friend StandardOutputStream& operator<<(StandardOutputStream&, std::_Smanip<int>);
#endif
  
  // Overloaded output functions
  StandardOutputStream& operator<<(const std::string& output);
  StandardOutputStream& operator<<(const char* output);
  StandardOutputStream& operator<<(char output);
  StandardOutputStream& operator<<(int output);
  StandardOutputStream& operator<<(unsigned int output);
  StandardOutputStream& operator<<(short int output);
  StandardOutputStream& operator<<(unsigned short int output);
  StandardOutputStream& operator<<(long int output);
  StandardOutputStream& operator<<(unsigned long int output);
  StandardOutputStream& operator<<(long long int output);
  StandardOutputStream& operator<<(float output);
  StandardOutputStream& operator<<(double output);
  StandardOutputStream& operator<<(long double output);
  StandardOutputStream& operator<<(bool output);


private:
  //! Hide copy constructor and =
  StandardOutputStream(const StandardOutputStream&) {}
  void operator=(const StandardOutputStream&) {}

protected:
  // The universal data-writer. All the write functions call this
  template<typename T>
  StandardOutputStream& writePrimitive(T output);

  // Get the internal ostream
  std::ostream& getOstream() {return *os;}

  //! Is the stream open?
  bool is_open() {return open;}

private:
  std::ostream* os;
  bool open;
};

// NOTE: A write of an array is *NOT* provided since these are usual special purpose
//template<class T>
//inline
//StandardOutputStream& operator<<(StandardOutputStream& s, const multi1d<T>& d);

// Make global (parallel) versions of stdin, stdout, stderr
// Put this in a different namespace to avoid collisions
namespace CCIO{
  extern StandardInputStream  cin;
  extern StandardOutputStream cout;
  extern StandardOutputStream cerr;
  void header(std::string name);
}



#endif //COMM_IO_H_

