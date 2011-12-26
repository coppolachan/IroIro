/*! 
 * @file comm_io.cpp 
 * @brief Definition of MPI safe input/output routines
 *
 * These are the parallel versions of the STL cout,cerr,cin.
 * The user can control which nodes are involved in output/input.
 */


#include "comm_io.hpp"
#include "communicator.h"


/*!
 *@brief CommonCode Input/Output overloaded functions
 *
 *Standard-IO-like input and output operations 
 *
 * Credits to the Chroma team
 */
namespace CCIO
{
  StandardInputStream  cin;
  StandardOutputStream cout;
  StandardOutputStream cerr;

  void header(std::string name){
    cout << "\n--------------------------------\n";
    cout << name;
    cout << "\n--------------------------------\n";

  }


}


//-----------------------------------------
//! stdin support

StandardInputStream::StandardInputStream() {open=false; is=0;}

StandardInputStream::~StandardInputStream() {}

void StandardInputStream::init(std::istream* b)
{
  is = b;
  open = true;
}
/*
// Propagate status to all nodes
bool StandardInputStream::fail()
{
  bool s;

  if (Communicator::instance()->primaryNode()) 
    s = getIstream().fail();

  Internal::broadcast(s);
  return s;
}

// String reader
StandardInputStream& StandardInputStream::operator>>(std::string& input)
{
  // Only primary node can grab string
  if (Communicator::instance()->primaryNode()) 
  {
    getIstream() >> input;
  }

  // broadcast string
  Internal::broadcast_str(input);

  return *this;
}

// Readers
StandardInputStream& StandardInputStream::operator>>(char& input) 
{
  return readPrimitive<char>(input);
}
StandardInputStream& StandardInputStream::operator>>(int& input) 
{
  return readPrimitive<int>(input);
}
StandardInputStream& StandardInputStream::operator>>(unsigned int& input)
{
  return readPrimitive<unsigned int>(input);
}
StandardInputStream& StandardInputStream::operator>>(short int& input)
{
  return readPrimitive<short int>(input);
}
StandardInputStream& StandardInputStream::operator>>(unsigned short int& input)
{
  return readPrimitive<unsigned short int>(input);
}
StandardInputStream& StandardInputStream::operator>>(long int& input)
{
  return readPrimitive<long int>(input);
}
StandardInputStream& StandardInputStream::operator>>(unsigned long int& input)
{
  return readPrimitive<unsigned long int>(input);
}
StandardInputStream& StandardInputStream::operator>>(long long int& input)
{
  return readPrimitive<long long int>(input);
}
StandardInputStream& StandardInputStream::operator>>(float& input)
{
  return readPrimitive<float>(input);
}
StandardInputStream& StandardInputStream::operator>>(double& input)
{
  return readPrimitive<double>(input);
}
StandardInputStream& StandardInputStream::operator>>(long double& input)
{
  return readPrimitive<long double>(input);
}
StandardInputStream& StandardInputStream::operator>>(bool& input)
{
  return readPrimitive<bool>(input);
}

template< typename T>
StandardInputStream& StandardInputStream::readPrimitive(T& input)
{
  if (Communicator::instance()->primaryNode())
    getIstream() >> input;

  // Now broadcast back out to all nodes
  Internal::broadcast(input);

  return *this;
}
*/

//-----------------------------------------
//! stdout support
StandardOutputStream::StandardOutputStream() {open=false; os=0;}

void StandardOutputStream::init(std::ostream* b)
{
  os = b;
  open = true;
}

void StandardOutputStream::flush()
{
  if (is_open()) 
  {
    if (Communicator::instance()->primaryNode()) 
      getOstream().flush();
  }
}

// Propagate status to all nodes
/*
bool StandardOutputStream::fail()
{
  bool s;

  if (Communicator::instance()->primaryNode()) 
    s = getOstream().fail();

  Internal::broadcast(s);
  return s;
}
*/

StandardOutputStream::~StandardOutputStream() {}


// Hook for manipulators
StandardOutputStream& StandardOutputStream::operator<<(std::ostream& (*op)(std::ostream&))
{
  // Call the function passed as parameter with this stream as the argument.
  // Hmm, for now do this on only the primary stream. This could potentially
  // cause problems since somebody might set hex mode or something and it
  // won't be done on all nodes.
  if (Communicator::instance()->primaryNode())
    (*op)(getOstream());

  return *this;
}


  

  StandardOutputStream& operator<<(StandardOutputStream& obj, std::_Setw f) {
    if (Communicator::instance()->primaryNode()) 
      obj.getOstream() << f;
    
    return obj;
  }
  
  StandardOutputStream& operator<<(StandardOutputStream& obj, std::_Setprecision f) {
    if (Communicator::instance()->primaryNode()) 
      obj.getOstream() << f;
    
    return obj;
  }


StandardOutputStream& StandardOutputStream::operator<<(const std::string& output)
{

  if (Communicator::instance()->primaryNode())
    getOstream() << output;

  return *this;
}

StandardOutputStream& StandardOutputStream::operator<<(const char* output)
{
  if (Communicator::instance()->primaryNode())
    getOstream() << output;

  return *this;
}

StandardOutputStream& StandardOutputStream::operator<<(char output) 
{
  return writePrimitive<char>(output);
}

StandardOutputStream& StandardOutputStream::operator<<(int output) 
{
  return writePrimitive<int>(output);
}

StandardOutputStream& StandardOutputStream::operator<<(unsigned int output)
{
  return writePrimitive<unsigned int>(output);
}

StandardOutputStream& StandardOutputStream::operator<<(short int output)
{
  return writePrimitive<short int>(output);
}

StandardOutputStream& StandardOutputStream::operator<<(unsigned short int output)
{
  return writePrimitive<unsigned short int>(output);
}

StandardOutputStream& StandardOutputStream::operator<<(long int output)
{
  return writePrimitive<long int>(output);
}

StandardOutputStream& StandardOutputStream::operator<<(unsigned long int output)
{
  return writePrimitive<unsigned long int>(output);
}

StandardOutputStream& StandardOutputStream::operator<<(long long int output)
{
  return writePrimitive<long long int>(output);
}

StandardOutputStream& StandardOutputStream::operator<<(float output)
{
  if (Communicator::instance()->primaryNode())
  {
    std::streamsize initPrec = getOstream().precision();
    getOstream().precision(7);
    getOstream() << output;
    getOstream().precision(initPrec);
  }

  return *this;
}

StandardOutputStream& StandardOutputStream::operator<<(double output)
{
  if (Communicator::instance()->primaryNode())
  {
    std::streamsize initPrec = getOstream().precision();
    getOstream().precision(15);
    getOstream() << output;
    getOstream().precision(initPrec);
  }

  return *this;
}

StandardOutputStream& StandardOutputStream::operator<<(long double output)
{
  return writePrimitive<long double>(output);
}


StandardOutputStream& StandardOutputStream::operator<<(bool output)
{
  return writePrimitive<bool>(output);
}

template<typename T>
StandardOutputStream& StandardOutputStream::writePrimitive(T output)
{
  if (Communicator::instance()->primaryNode())
    getOstream() << output;

  return *this;
}

